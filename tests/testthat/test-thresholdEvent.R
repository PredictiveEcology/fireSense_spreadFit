## The `estimateThreshold` and `debug` events: what they hand to runSpreadWithoutDEoptim(), which is
## replaced in the module's environment by a recorder with the same arguments.

recordingRSWD <- function(rec, value) {
  f <- function() {
    rec$args <- as.list(environment())
    value
  }
  formals(f) <- formals(runSpreadWithoutDEoptim)
  f
}

test_that("without SNLL_FS_thresh the threshold is calibrated, seeded by the polygon id", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE, SNLL_FS_thresh = NULL, iterThresh = 6L, objfunFireReps = 8L))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, 777))
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_identical(rec$args$seed, .elfSeed("9.9"))
  expect_identical(rec$args$iterThresh, 6L)
  expect_identical(rec$args$objfunFireReps, 8L)
  expect_identical(rec$args$tests, c("adTest", "SNLL_FS"))
  expect_identical(rec$args$mode, "fit")
  expect_identical(rec$args$lower, SpaDES.core::params(sim)[[moduleName]]$lower)
  expect_identical(rec$args$annualDTx1000$year2002$CMDsm, c(12000L, 22000L, 32000L))
  ## and the calibrated value is the threshold of the fit
  expect_identical(as.vector(fitRec$deArgs$thresh), 777)   # Cache() adds attributes
})

test_that("a supplied SNLL_FS_thresh is used as is: no calibration", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE, SNLL_FS_thresh = 250L))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, 777))
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_null(rec$args)
  expect_identical(fitRec$deArgs$thresh, 250L)
})

test_that("debug mode evaluates the objective function instead of fitting", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE, mode = "debug", iterThresh = 3L, DEoptimTests = c("SNLL_FS", "adTest")))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, NULL))
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_null(fitRec$deArgs)                             # no DEoptim
  expect_null(fitRec$geoArgs)                            # nothing written to the ledger
  expect_identical(rec$args$mode, "debug")
  expect_identical(rec$args$iterThresh, 3L)
  expect_identical(rec$args$tests, c("SNLL_FS", "adTest"))
  expect_null(rec$args$seed)                             # the debug event does not fix the seed
  expect_identical(rec$args$covMinMax, sim$covMinMax_spread)
  done <- SpaDES.core::completed(sim)
  expect_identical(tail(done$eventType, 1), "debug")
})
