## DEoptim settings reach DEoptim. Eliot, 2026-09-15: "Any user passed args should pass into the
## DEoptim processes." The module passed strategy, trace, initialpop and .c only, and .c did not
## reach DEoptim at all; `DEoptimControl` carries any other DEoptim.control() setting (CR, F, p,
## reltol, ...) through fireSenseUtils::runDEoptim() to clusters::clusterSetup().

test_that("DEoptimControl is a list parameter; by default it sets only p, for strategy 6", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  p <- md$parameters
  i <- which(p$paramName == "DEoptimControl")
  expect_length(i, 1L)
  expect_identical(p$paramClass[[i]], "list")
  expect_identical(p$default[[i]], list(p = 0.1))
})

test_that("the fit passes DEoptimControl, strategy and .c to runDEoptim", {
  sim <- toySim(list(DEoptimControl = list(p = 0.2, CR = 0.7), strategy = 3L, .c = 0.1))
  rec <- mockFitAndLedger(sim)
  suppressMessages(SpaDES.core::spades(sim))
  ## the parameter as given (a starting population is added only with passing sets or on the retry after a
  ## `deoptimNoFeasibleMember` stop, test-seedPopulation.R; toySim() has none)
  expect_identical(rec$deArgs$DEoptimControl, list(p = 0.2, CR = 0.7))
  expect_identical(rec$deArgs$strategy, 3L)
  expect_identical(rec$deArgs$.c, 0.1)
})
