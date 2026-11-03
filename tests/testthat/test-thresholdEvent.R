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

test_that("without firstBlockThresh the threshold is calibrated, seeded by the polygon id", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(firstBlockThresh = NULL, iterThresh = 6L, objfunFireReps = 8L))
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

test_that("a supplied firstBlockThresh is used as is: no calibration", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(firstBlockThresh = 250))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, 777))
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_null(rec$args)
  expect_identical(fitRec$deArgs$thresh, 250)
})

test_that("debug mode evaluates the objective function instead of fitting", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(mode = "debug", iterThresh = 3L, DEoptimTests = c("SNLL_FS", "adTest")))
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

test_that("one option, fireSenseUtils.assertions, sets the objective's diagnostics for fit and threshold", {
  ## off by default (they are slow); the threshold calibration takes the module parameter, not its own default
  expect_identical(formals(runSpreadWithoutDEoptim)$doObjFunAssertions, quote(expr = ))   # no default of its own
  for (opt in c(FALSE, TRUE)) {
    withr::local_options(fireSenseUtils.assertions = opt)
    rec <- new.env(); fitRec <- new.env()
    sim <- toySim(list(firstBlockThresh = NULL, iterThresh = 3L, objfunFireReps = 2L))
    expect_identical(SpaDES.core::P(sim, "doObjFunAssertions", moduleName), opt)
    mockFitAndLedger(sim, fitRec)
    mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, 777))
    suppressMessages(SpaDES.core::spades(sim))
    expect_identical(rec$args$doObjFunAssertions, opt)
  }
})

test_that("an objective without an early stop skips the threshold calibration", {
  fitRec <- new.env(); calls <- 0L
  sim <- toySim(list(firstBlockThresh = NULL, DEoptimTests = "crps"))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = function(...) { calls <<- calls + 1L; stop("the calibration must not run") })
  expect_message(sim <- SpaDES.core::spades(sim), "no early stop")
  expect_identical(calls, 0L)
  expect_identical(fitRec$deArgs$thresh, Inf)
  expect_null(fitRec$deArgs$initialpop)
})

test_that("an objective with an early stop is still calibrated, once", {
  rec <- new.env(); fitRec <- new.env(); calls <- 0L
  sim <- toySim(list(firstBlockThresh = NULL, DEoptimTests = c("adTest", "SNLL_FS")))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = function(...) { calls <<- calls + 1L; 777 })
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_identical(calls, 1L)
  expect_identical(as.vector(fitRec$deArgs$thresh), 777)
})
