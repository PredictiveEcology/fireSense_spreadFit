## The whole fit path -- init, spreadFitPrepare, estimateThreshold, run -- on toy inputs, with the two
## calls that leave the machine replaced by recorders (see helper-toyInputs.R): runDEoptim() returns
## a small fixed DEoptim result, CacheGeo() returns the row it was asked to write.

m <- function(sim) SpaDES.core::params(sim)[[moduleName]]

fitted <- function(params = list(), objects = list(), rec = new.env()) {
  sim <- toySim(c(list(stopIfNoPreRunFit = FALSE), params), objects)
  mockFitAndLedger(sim, rec)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  list(sim = sim, rec = rec)
}

test_that("every scheduled event runs, in order, and nothing is left queued", {
  out <- fitted()
  done <- SpaDES.core::completed(out$sim)
  expect_identical(done$eventType[done$moduleName == moduleName],
                   c(".inputObjects", "init", "spreadFitPrepare", "spreadFitPrepare",
                     "estimateThreshold", "run", "postFitDiagnostics"))
  expect_length(queued(out$sim), 0L)
})

test_that("runDEoptim() receives the covariates as integers x 1000, in parameter order", {
  a <- fitted()$rec$deArgs
  expect_identical(names(a$annualDTx1000), c("year2001", "year2002"))
  expect_equal(a$annualDTx1000$year2001,
               data.frame(pixelID = c(2L, 3L, 4L, 55L, 56L, 57L),
                          CMDsm = c(10000L, 20000L, 30000L, 40000L, 25000L, 15000L),
                          youngAge = c(0L, 1000L, 0L, 0L, 0L, 1000L)))
  expect_identical(class(a$annualDTx1000$year2001), "data.frame")   # not a data.table
  expect_type(a$annualDTx1000$year2001$CMDsm, "integer")
  ## fuel reaches the fit as LINEAR biomass x 1000, not as the log it was supplied on: 0 on the
  ## logMinB() floor is exactly 0, and 22601 is 22601000
  expect_equal(a$nonAnnualDTx1000$year2001_year2002$class1,
               c(0, 300, 1000, 5000, 7803, 22601, 300, 1000, 5000) * 1000, tolerance = 1e-9)
  expect_equal(a$nonAnnualDTx1000$year2001_year2002$class2,
               c(40, 0, 2000, 0, 12000, 0, 500, 0, 9000) * 1000, tolerance = 1e-9)
  expect_type(a$nonAnnualDTx1000$year2001_year2002$class1, "integer")
  expect_identical(a$nonAnnualDTx1000$year2001_year2002$nf,                # not fuel: unchanged
                   c(0L, 100L, 200L, 300L, 400L, 500L, 600L, 700L, 800L))
  expect_identical(a$nonAnnualDTx1000$year2001_year2002$pixelID,           # pixelID is not scaled
                   c(2L, 3L, 4L, 55L, 56L, 57L, 88L, 89L, 90L))
  expect_equal(a$fireBufferedListDT$year2002,
               data.frame(pixelID = c(88L, 89L, 90L), buffer = c(0L, 1L, 0L), ids = 21L))
  expect_equal(a$historicalFires$year2001,
               data.frame(size = c(5L, 12L), date = "year2001", ids = 11:12, cells = c(3L, 56L)))
})

test_that("the simList's own covariate tables are not multiplied by 1000", {
  sim <- fitted()$sim
  expect_equal(sim$fireSense_annualSpreadFitCovariates$year2001$CMDsm, c(10, 20, 30, 40, 25, 15))
  expect_equal(sim$fireSense_nonAnnualSpreadFitCovariates[[1]]$nf, seq(0, 0.8, by = 0.1))
})

test_that("runDEoptim() receives the parameters, bounds, threshold and formula", {
  out <- fitted(list(iterDEoptim = 40L, nCoresNeeded = 12L, cores = c("hostA", "hostB"),
                     objfunFireReps = 9L, .rep = 3L, DEoptimControl = list(CR = 0.7), .c = 0.4,
                     SNLL_FS_thresh = 321L))
  a <- out$rec$deArgs
  expect_identical(a$itermax, 40L)
  expect_identical(a$nCoresNeeded, 12L)
  expect_identical(a$cores, c("hostA", "hostB"))
  expect_identical(a$Nreps, 9L)
  expect_identical(a$rep, 3L)
  expect_identical(a$DEoptimControl, list(CR = 0.7))
  expect_identical(a$.c, 0.4)
  expect_identical(a$thresh, 321L)                       # SNLL_FS_thresh, not an estimate
  expect_identical(a$runName, "toyRun")
  expect_identical(a$formulaToFit, "~ 0 + CMDsm + youngAge + class1 + class2 + nf")
  expect_identical(a$lower, m(out$sim)$lower)
  expect_identical(names(a$upper), c("maxAsymptote", "inflectionPoint1",
                                     "CMDsm", "youngAge", "class1", "class2", "nf", "yearSpreadSD"))
  expect_identical(a$mutuallyExclusive, list(youngAge = c("class", "nonForest", "class1", "class2", "nf")))
  expect_identical(a$covMinMax, out$sim$covMinMax_spread)
  expect_identical(a$tests, c("adTest", "SNLL_FS"))
  expect_identical(a$maxFireSpread, 0.28)
})

test_that("runDEoptim() always receives iterStep = 1, even if a user still sets the (now unknown) parameter", {
  a <- suppressWarnings(fitted(list(iterStep = 10L)))$rec$deArgs
  expect_identical(a$iterStep, 1L)
})

## 2026-09-28: DEoptim progress figures drawn after every generation took 16% of a fit's wall time;
## .plotInterval (default 25) sets how many generations pass between them.
test_that("runDEoptim() receives .plotInterval as plotEvery", {
  expect_identical(fitted()$rec$deArgs$plotEvery, 25)
  expect_identical(fitted(list(.plotInterval = 5))$rec$deArgs$plotEvery, 5)
})

test_that("years without fires, or without a fire buffer, are dropped from every annual list", {
  objs <- toyObjects()
  dt <- data.table::data.table
  ## year2003: covariates only (no buffer, no fire). year2004: covariates and buffer, but no fire.
  objs$fireSense_annualSpreadFitCovariates$year2003 <- dt(pixelID = 1:2, CMDsm = c(11, 12), youngAge = c(0, 0))
  objs$fireSense_annualSpreadFitCovariates$year2004 <- dt(pixelID = 1:2, CMDsm = c(13, 14), youngAge = c(0, 0))
  objs$fireBufferedListDT$year2004 <- dt(pixelID = 1:2, buffer = c(1L, 0L), ids = 41L)
  a <- fitted(objects = objs[c("fireSense_annualSpreadFitCovariates", "fireBufferedListDT")])$rec$deArgs
  expect_identical(names(a$annualDTx1000), c("year2001", "year2002"))
  expect_identical(names(a$fireBufferedListDT), c("year2001", "year2002"))
  expect_identical(names(a$historicalFires), c("year2001", "year2002"))
})

test_that("sim$DE is the DEoptim result ordered by each block's best value", {
  sim <- fitted()$sim
  ## toyDE() blocks have bestvalit 30, 10, 20
  expect_identical(vapply(sim$DE, function(d) d$member$bestvalit, numeric(1)), c(10, 20, 30))
})

test_that("the ledger row: polygon id, 5 best members, their values, and the covariate ranges", {
  out <- fitted()
  row <- out$sim$studyAreaWithSpreadParams
  expect_identical(nrow(row), 1L)
  expect_identical(names(row), c("geometry", "numIterations", "objFunVal", "params", "sppEquiv",
                                 "nonForestedLCCGroups", "missingLCCgroup", "covMinMax_spread",
                                 "polygonID"))
  expect_identical(row$polygonID, "9.9")                 # .ELFind, not .runName
  expect_identical(row$numIterations[[1]], 3L)           # toyDE() has 3 blocks
  ## toyDE(): member k has first parameter k and value 8 - k, so the 5 best are members 7..3
  best <- row$params[[1]]
  ## hillSlope1 is fixed at 1, not fitted, so it is not in names(m(out$sim)$lower); the ledger row
  ## gets it back, right after maxAsymptote, so it predicts like an old row (addHillSlope1ToLedger())
  expect_identical(names(best), append(names(m(out$sim)$lower), "hillSlope1", after = 1L))
  expect_equal(best$maxAsymptote, c(7, 6, 5, 4, 3))
  expect_identical(best$hillSlope1, rep(1, 5))
  expect_equal(row$objFunVal[[1]], c(1, 2, 3, 4, 5))
  expect_identical(row$covMinMax_spread[[1]], out$sim$covMinMax_spread)
  expect_true(sf::st_equals(row, toyStudyArea(), sparse = FALSE)[1, 1])
  expect_identical(sf::st_crs(row), sf::st_crs(toyCRS))
})

test_that("the ledger is written to the configured Drive folder and file, as an update", {
  g <- fitted(list(spreadFitGoogleDriveFolder = "https://drive.google.com/drive/folders/someFolder",
                   spreadFitFilename = "toyLedger.rds"))$rec$geoArgs
  expect_identical(g$cloudFolderID, "https://drive.google.com/drive/folders/someFolder")
  expect_identical(g$targetFile, "toyLedger.rds")
  expect_identical(g$action, "update")
  expect_identical(g$destinationPath, toyPaths()$inputPath)
  expect_identical(g$studyAreaFireSense$polygonID, "9.9")
  expect_true(sf::st_equals(g$domain, toyStudyArea(), sparse = FALSE)[1, 1])
})

test_that("by default (\"latest\") the ledger file is named for the fit's fire years and model", {
  g <- fitted()$rec$geoArgs
  ## the model tag (e.g. "_linearFuel_esc50") is fireSenseUtils' to set; only the years come from here
  expect_identical(g$targetFile, paste0("fireSenseParams_2001-2002", fireSenseUtils::spreadFitFileTag, ".rds"))  # the toy's covariate years
  expect_identical(ledgerWriteFile("latest", 1985:2024, c("year2001", "year2002")),
                   paste0("fireSenseParams_1985-2024", fireSenseUtils::spreadFitFileTag, ".rds"))  # dataPrepFit's window wins
  expect_identical(ledgerWriteFile("mine.rds", 1985:2024, "year2001"), "mine.rds")
})

test_that("sppEquiv and the land-cover groups in the simList are recorded in the row", {
  out <- fitted(objects = list(sppEquiv = data.frame(sp = c("a", "b")),
                               nonForestedLCCGroups = list(nf = 1:3), missingLCCgroup = "nf"))
  row <- out$sim$studyAreaWithSpreadParams
  expect_identical(row$sppEquiv[[1]], data.frame(sp = c("a", "b")))
  expect_identical(row$nonForestedLCCGroups[[1]], list(nf = 1:3))
  expect_identical(row$missingLCCgroup[[1]], "nf")
})

test_that("an .ELFind that cannot key the ledger stops before anything is written", {
  rec <- new.env()
  expect_error(fitted(objects = list(.ELFind = c("9.9", "9.8")), rec = rec),
               "`sim\\$.ELFind` must be a single non-empty character identifying the polygon being fit; got: 9.9, 9.8")
  expect_null(rec$geoArgs)
  expect_error(fitted(objects = list(.ELFind = "")), "must be a single non-empty character")
})

test_that("stale spreadFitAdditionalColNames are replaced by fireSenseUtils' names", {
  sim <- fitted(objects = list(spreadFitAdditionalColNames = c("a", "b")))$sim
  expect_identical(sim$spreadFitAdditionalColNames, fireSenseUtils::spreadFitAdditionalColNamesTxt)
})

test_that("visualizeDEoptim is pointed at the module's figure folder", {
  out <- fitted(list(visualizeDEoptim = file.path(tempdir(), "elsewhere")))
  expect_identical(basename(out$rec$deArgs$visualizeDEoptim), moduleName)
  expect_identical(basename(m(out$sim)$visualizeDEoptim), moduleName)
})

test_that("with a ledger row for this polygon nothing is fitted or written", {
  rec <- new.env()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE), list(studyAreaWithSpreadParams = toyLedger("9.9")))
  mockFitAndLedger(sim, rec)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_null(rec$deArgs)
  expect_null(rec$geoArgs)
  expect_null(sim$DE)
  expect_identical(sim$studyAreaWithSpreadParams$polygonID, "9.9")
})

test_that("refitExisting fits and writes although the ledger has this polygon", {
  rec <- new.env()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE, refitExisting = TRUE),
                list(studyAreaWithSpreadParams = toyLedger("9.9")))
  mockFitAndLedger(sim, rec)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_identical(rec$deArgs$itermax, 5000L)
  expect_identical(rec$geoArgs$studyAreaFireSense$polygonID, "9.9")
  expect_equal(sim$studyAreaWithSpreadParams$params[[1]]$maxAsymptote, c(7, 6, 5, 4, 3))
})
