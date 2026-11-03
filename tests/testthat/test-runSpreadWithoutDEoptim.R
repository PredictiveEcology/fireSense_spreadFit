## runSpreadWithoutDEoptim(): the `debug` event's evaluation of the objective, which is replaced by a
## recorder (the real one simulates fires and belongs to fireSenseUtils).

fires <- list(year2001 = data.frame(size = c(5L, 12L), date = "year2001", ids = 11:12, cells = c(3L, 56L)),
              year2002 = data.frame(size = c(3L, 1L), date = "year2002", ids = 21:22, cells = c(89L, 12L)))
lower <- c(a = 0, b = 10)
upper <- c(a = 1, b = 20)

callIt <- function(..., nSets = 4L, seed = 42L, lowerBounds = lower, upperBounds = upper) {
  utils::capture.output(out <- runSpreadWithoutDEoptim(
    nSets = nSets, lower = lowerBounds, upper = upperBounds, fireSense_spreadFormula = "~ 0 + a + b",
    flammableRTM = "theLandscape", annualDTx1000 = list(), nonAnnualDTx1000 = list(),
    fireBufferedListDT = list(), historicalFires = fires, covMinMax = NULL, objfunFireReps = 7L,
    maxFireSpread = 0.28, tests = "crps", formulaToFit = "~ 0 + a + b", seed = seed,
    doObjFunAssertions = FALSE, ...))
  out
}

## records each call
recorder <- function(rec) {
  function(par, ...) {
    rec$calls[[length(rec$calls) + 1L]] <- c(list(par = par), list(...))
    1000
  }
}

test_that("evaluates every random parameter set in turn and returns NULL", {
  rec <- new.env()
  local_mocked_bindings(.objfunSpreadFit = recorder(rec))
  expect_null(callIt())
  expect_length(rec$calls, 4L)
  pars <- do.call(rbind, lapply(rec$calls, `[[`, "par"))
  expect_true(all(pars[, 1] >= 0 & pars[, 1] <= 1))      # each parameter inside its own bounds
  expect_true(all(pars[, 2] >= 10 & pars[, 2] <= 20))
  expect_false("thresh" %in% names(rec$calls[[1]]))         # no early stop
  ## what the objective function is told
  one <- rec$calls[[1]]
  expect_identical(one$landscape, "theLandscape")
  expect_identical(one$Nreps, 7L)
  expect_identical(one$tests, "crps")
  expect_identical(one$maxFireSpread, 0.28)
  expect_identical(one$historicalFires, fires)
  expect_true(one$weighted)
  expect_true(one$plot.it)
})

test_that("the seed fixes the draws: parameter draws, R >= 3.6 sampling", {
  rec <- new.env()
  local_mocked_bindings(.objfunSpreadFit = recorder(rec))
  callIt(seed = 42L)
  first <- rec$calls
  rec$calls <- NULL
  callIt(seed = 42L)
  expect_identical(rec$calls, first)
  expect_equal(first[[1]]$par, c(a = 0.914806043496355, b = 19.3707541329786))
  rec$calls <- NULL
  callIt(seed = 43L)
  expect_false(identical(rec$calls[[1]]$par, first[[1]]$par))
})

test_that("supplied `pars` are evaluated as given", {
  rec <- new.env()
  local_mocked_bindings(.objfunSpreadFit = recorder(rec))
  callIt(pars = c(0.5, 15))
  expect_length(rec$calls, 1L)
  ## named with `names(lower)`: unnamed and the same length, so the objective function can tell
  ## apart a trailing yearSpreadSD from a logistic parameter (see the naming test below)
  expect_identical(rec$calls[[1]]$par, c(a = 0.5, b = 15))
})

test_that("drawn parameter sets are named, so a trailing yearSpreadSD bound is recognised", {
  ## since yearSpreadSD became a default trailing bound, fireSenseUtils:::.objfunSpreadFit tells it apart from a
  ## logistic parameter only via names(par): an unnamed draw made every evaluation error with
  ## "logistic with 4 parameters not tested yet"
  rec <- new.env()
  local_mocked_bindings(.objfunSpreadFit = recorder(rec))
  callIt(lowerBounds = c(lower, yearSpreadSD = 0), upperBounds = c(upper, yearSpreadSD = 1))
  expect_identical(names(rec$calls[[1]]$par), c("a", "b", "yearSpreadSD"))
})
