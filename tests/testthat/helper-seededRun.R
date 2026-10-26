## The run event with the threshold calibration's passing sets, for the tests of DEoptim's starting population
## (test-seedPopulation.R, test-emulator.R).

## three sets inside the bounds, columns as DEoptim uses them (the bounds are made by spreadFitPrepare, so the
## mocked calibration builds them from the bounds it is given)
passingSets <- function(lo, up) {
  m <- rbind(lo + 0.25 * (up - lo), lo + 0.5 * (up - lo), lo + 0.75 * (up - lo))
  colnames(m) <- names(lo)
  m
}

noFeasible <- function() structure(class = c("deoptimNoFeasibleMember", "error", "condition"),
                                   list(message = "every member of the population holds the fail value", call = NULL))

## The run event with the calibration and runDEoptim() replaced. `calibration` is what the cached
## calibration returns; `failFirst` makes the first fit stop as clusters does; `mocks`, a function of the record
## `rec`, returns further functions of the module to replace (mockInModule()).
seededRun <- function(calibration, failFirst = FALSE, params = list(), mocks = function(rec) list()) {
  ## the module's Cache() of the calibration would serve the previous test's mocked result
  withr::local_options(reproducible.useCache = FALSE)
  rec <- new.env(); rec$calls <- list()
  sim <- toySim(c(list(firstBlockThresh = NULL), params))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runSpreadWithoutDEoptim = function(lower, upper, ...) {
      rec$lower <- lower; rec$upper <- upper; rec$passing <- passingSets(lower, upper)
      if (is.function(calibration)) calibration(rec) else calibration
    },
    runDEoptim = function(...) {
      args <- list(...)
      rec$calls <- c(rec$calls, list(args))
      if (failFirst && length(rec$calls) == 1L) stop(noFeasible())
      toyDE(length(args$lower))
    },
    feasibleInitialPopulation = function(n, lower, upper, gates, seed, ...) {
      m <- matrix(rep(lower, each = n) + 0.1, n, length(lower), dimnames = list(NULL, names(lower)))
      attr(m, "feasible") <- 5L
      m
    })
  do.call(mockInModule, c(list(sim), mocks(rec)))
  suppressMessages(SpaDES.core::spades(sim))
  rec
}
