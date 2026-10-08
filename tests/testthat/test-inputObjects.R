## .inputObjects(): the defaults that need no download. The `studyArea` and `rasterToMatch` defaults
## download and are not tested (and `rasterToMatch` has no sourceURL to download from).

test_that("fireSense_spreadFormula is mandatory", {
  expect_error(toySim(objects = list(fireSense_spreadFormula = NULL)),
               "fireSense_spreadFormula must be supplied")
})

test_that(".ELFind falls back to .runName, and says so", {
  msgs <- testthat::capture_messages(
    sim <- SpaDES.core::simInit(times = list(start = 0, end = 1), modules = moduleName, paths = toyPaths(),
                                params = stats::setNames(list(list(.useCache = FALSE)), moduleName),
                                objects = within(toyObjects(), rm(.ELFind))))
  expect_identical(sim$.ELFind, "toyRun")
  expect_match(paste(msgs, collapse = ""), "`.ELFind` not supplied; keying the shared fit ledger on `.runName` \\('toyRun'\\)")
})

test_that("a supplied .ELFind is kept", {
  sim <- toySim()
  expect_identical(sim$.ELFind, "9.9")
  expect_identical(sim$.runName, "toyRun")
})

test_that("spreadFitAdditionalColNames defaults to the ledger's list-columns", {
  sim <- toySim()
  expect_identical(sim$spreadFitAdditionalColNames,
                   c("numIterations", "objFunVal", "params", "sppEquiv", "nonForestedLCCGroups",
                     "missingLCCgroup", "covMinMax_spread"))
  sim <- toySim(objects = list(spreadFitAdditionalColNames = c("a", "b")))
  expect_identical(sim$spreadFitAdditionalColNames, c("a", "b"))   # not overwritten here (run resets it)
})

test_that("supplied inputs are not replaced", {
  sim <- toySim()
  expect_identical(terra::ncell(sim$rasterToMatch), 100)
  expect_equal(as.numeric(sf::st_area(sim$studyArea)), 1e6)        # 1000 m x 1000 m
})
