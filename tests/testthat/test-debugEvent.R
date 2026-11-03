## The `debug` event: what it hands to runSpreadWithoutDEoptim(), which is replaced in the module's
## environment by a recorder with the same arguments.

test_that("debug mode evaluates the objective function instead of fitting", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(mode = "debug", objfunFireReps = 8L, DEoptimTests = c("adTest", "pitAD")))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, NULL))
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_null(fitRec$deArgs)                             # no DEoptim
  expect_null(fitRec$geoArgs)                            # nothing written to the ledger
  expect_identical(rec$args$objfunFireReps, 8L)
  expect_identical(rec$args$tests, c("adTest", "pitAD"))
  expect_identical(rec$args$lower, SpaDES.core::params(sim)[[moduleName]]$lower)
  expect_identical(rec$args$annualDTx1000$year2002$CMDsm, c(12000L, 22000L, 32000L))
  expect_null(rec$args$seed)                             # the debug event does not fix the seed
  expect_identical(rec$args$covMinMax, sim$covMinMax_spread)
  done <- SpaDES.core::completed(sim)
  expect_identical(tail(done$eventType, 1), "debug")
})

test_that("one option, fireSenseUtils.assertions, sets the objective's diagnostics for the fit and the debug call", {
  ## off by default (they are slow); runSpreadWithoutDEoptim() takes the module parameter, not its own default
  expect_identical(formals(runSpreadWithoutDEoptim)$doObjFunAssertions, quote(expr = ))   # no default of its own
  for (opt in c(FALSE, TRUE)) {
    withr::local_options(fireSenseUtils.assertions = opt)
    expect_identical(SpaDES.core::P(toySim(), "doObjFunAssertions", moduleName), opt)
    expect_identical(debugArgs()$doObjFunAssertions, opt)
    expect_identical(fitArgs()$doObjFunAssertions, opt)
  }
})
