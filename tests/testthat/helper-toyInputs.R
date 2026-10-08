## Synthetic inputs for the event-level tests. Nothing here touches the network, Google Drive or a
## cluster: the two calls that would (`runDEoptim()` and `CacheGeo()`) are replaced in the module's
## environment by `mockInModule()`.
##
## The landscape is a 10 x 10 raster of 100 m (= 1 ha) pixels, so a fire of `POLY_HA = 5` is 5 pixels
## and the point (250, 950) is in row 1, column 3 = cell 3.

## pkgload sources helpers into the package namespace, where setup.R's `moduleName` and `testPaths`
## are not visible, so the helpers carry their own.
toyModule <- "fireSense_spreadFit"
toyEnv <- new.env()
toyPaths <- function() {
  if (is.null(toyEnv$paths)) {
    root <- withr::local_tempdir(.local_envir = testthat::teardown_env())
    ## tests run in tests/testthat: the directory holding the module directory is three levels up
    toyEnv$paths <- list(cachePath = file.path(root, "cache"), inputPath = file.path(root, "inputs"),
                         modulePath = normalizePath(file.path("..", "..", ".."), winslash = "/"),
                         outputPath = file.path(root, "outputs"))
    for (d in toyEnv$paths[c("cachePath", "inputPath", "outputPath")]) dir.create(d, recursive = TRUE)
  }
  toyEnv$paths
}

toyCRS <- "EPSG:3005"

toyRTM <- function() {
  terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000,
              crs = toyCRS, vals = 1L)
}

toyStudyArea <- function() {
  bb <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000), crs = sf::st_crs(toyCRS))
  sf::st_as_sf(sf::st_as_sfc(bb))
}

toyPoints <- function(x, y, id, yr, ha) {
  sf::st_as_sf(data.frame(x = x, y = y, FIRE_ID = id, YEAR = yr, POLY_HA = ha),
               coords = c("x", "y"), crs = toyCRS)
}

## A fresh copy every call: spreadFitPrep() changes some of its inputs by reference.
toyObjects <- function() {
  list(
    .runName = "toyRun", .ELFind = "9.9",
    rasterToMatch = toyRTM(), studyArea = toyStudyArea(),
    spreadFirePoints = list(
      year2001 = toyPoints(c(250, 550), c(950, 450), c(11L, 12L), 2001L, c(5, 12)), # cells 3 and 56
      year2002 = toyPoints(850, 150, 21L, 2002L, 3)),                               # cell 89
    fireBufferedListDT = list(
      year2001 = data.table::data.table(pixelID = c(2L, 3L, 4L, 55L, 56L, 57L),
                                        buffer = c(0L, 1L, 0L, 0L, 1L, 0L),
                                        ids = rep(c(11L, 12L), each = 3)),
      year2002 = data.table::data.table(pixelID = c(88L, 89L, 90L), buffer = c(0L, 1L, 0L), ids = 21L)),
    fireSense_annualSpreadFitCovariates = list(
      year2001 = data.table::data.table(pixelID = c(2L, 3L, 4L, 55L, 56L, 57L),
                                        CMDsm = c(10, 20, 30, 40, 25, 15),
                                        youngAge = c(0, 1, 0, 0, 0, 1)),
      year2002 = data.table::data.table(pixelID = c(88L, 89L, 90L),
                                        CMDsm = c(12, 22, 32), youngAge = c(0, 0, 1))),
    fireSense_nonAnnualSpreadFitCovariates = list(
      year2001_year2002 = data.table::data.table(
        pixelID = c(2L, 3L, 4L, 55L, 56L, 57L, 88L, 89L, 90L),
        ## fuel biomass as fireSense_dataPrepFit supplies it: fireSenseUtils::logMinB(), so 0 is on the floor
        class1 = fireSenseUtils::logMinB(c(0, 300, 1000, 5000, 7803, 22601, 300, 1000, 5000)),
        class2 = fireSenseUtils::logMinB(c(40, 0, 2000, 0, 12000, 0, 500, 0, 9000)),
        nf = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))), # cover-like: max <= 1
    fireSense_spreadFormula = "~ 0 + CMDsm + youngAge + class1 + class2 + nf"
  )
}

## A ledger as fireSense_dataPrepFit hands it over: one row per fitted polygon.
toyLedger <- function(ids) {
  sa <- toyStudyArea()
  sf::st_sf(polygonID = ids, geometry = rep(sf::st_geometry(sa), length(ids)))
}

## simInit() on the toy inputs. `params` and `objects` override the defaults; one set to NULL is
## removed (for a parameter, that means the module's default). `.useCache = FALSE`: the module caches `init` by default, and a cached `init` would
## hand one test the schedule another test made.
toySim <- function(params = list(), objects = list()) {
  ## not modifyList(): it merges data.tables column by column
  override <- function(base, new) {
    for (nm in names(new)) base[nm] <- list(new[[nm]])
    base[!vapply(base, is.null, logical(1))]
  }
  obj <- override(toyObjects(), objects)
  pars <- override(list(.useCache = FALSE, useCache_DE = FALSE, SNLL_FS_thresh = 500L), params)
  suppressMessages(
    SpaDES.core::simInit(times = list(start = 0, end = 1), modules = toyModule,
                         paths = toyPaths(), params = stats::setNames(list(pars), toyModule),
                         objects = obj))
}

## Run only the named events of this module.
runEvents <- function(sim, events) {
  suppressMessages(SpaDES.core::spades(sim, events = stats::setNames(list(events), toyModule)))
}

## Event types of this module still in the queue, in order.
queued <- function(sim) {
  ev <- SpaDES.core::events(sim)
  ev$eventType[ev$moduleName == toyModule]
}

## Replace a function the module calls, for this simList only.
mockInModule <- function(sim, ...) {
  fns <- list(...)
  for (nm in names(fns)) assign(nm, fns[[nm]], envir = sim$.mods[[toyModule]])
  invisible(sim)
}

## A DEoptim result as runDEoptim() returns it: one element per generation (iterStep is hard-coded to 1).
## `bestvalit` differs by generation so the ordering of sim$DE is checkable; the final population is 7 members whose first
## parameter is the member number, with values 7..1 so that members 7, 6, 5, 4, 3 are the 5 best.
toyDE <- function(nPar = 8L) {
  pop <- cbind(c(1, 2, 3, 4, 5, 6, 7), matrix(0, 7, nPar - 1L))
  gen <- function(bv) list(member = list(pop = pop, popval = c(7, 6, 5, 4, 3, 2, 1), bestvalit = bv))
  list(gen(30), gen(10), gen(20))
}

## The two calls of the `run` event that leave the machine, replaced by recorders.
mockFitAndLedger <- function(sim, rec = new.env()) {
  mockInModule(sim,
    runDEoptim = function(...) {
      rec$deArgs <- list(...)
      toyDE(length(rec$deArgs$lower))
    },
    ## `FUN` is never forced: in the module it is `le(studyAreaFireSense)`, which only evaluates
    ## inside the real CacheGeo(). The mock returns what the real one returns for a new ledger.
    CacheGeo = function(cloudFolderID, targetFile, domain, destinationPath, FUN, ..., action) {
      rec$geoArgs <- c(list(cloudFolderID = cloudFolderID, targetFile = targetFile, domain = domain,
                            destinationPath = destinationPath, action = action), list(...))
      list(...)$studyAreaFireSense
    })
  rec
}

## The module's source, parsed: the main file and its R/ helpers. For the tests that check a call by
## reading it, because running it needs a full simList.
moduleSource <- function() {
  files <- c(testthat::test_path("..", "..", "fireSense_spreadFit.R"),
             list.files(testthat::test_path("..", "..", "R"), "\\.R$", full.names = TRUE))
  do.call(c, lapply(files, parse, keep.source = FALSE))
}
