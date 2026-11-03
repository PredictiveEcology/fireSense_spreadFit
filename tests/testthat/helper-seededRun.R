## The run event, for the tests of DEoptim's starting population (test-seedPopulation.R, test-emulator.R).

## three sets inside the bounds, columns as DEoptim uses them (the bounds are made by spreadFitPrepare, so the
## tests build them from the bounds the fit is given)
passingSets <- function(lo, up) {
  m <- rbind(lo + 0.25 * (up - lo), lo + 0.5 * (up - lo), lo + 0.75 * (up - lo))
  colnames(m) <- names(lo)
  m
}

noFeasible <- function() structure(class = c("deoptimNoFeasibleMember", "error", "condition"),
                                   list(message = "every member of the population holds the fail value", call = NULL))

## The run event with runDEoptim() replaced; `rec` records its calls, and the bounds it was given. `failFirst`
## makes the first fit stop as clusters does; `mocks`, a function of the record `rec`, returns further functions
## of the module to replace (mockInModule()).
seededRun <- function(failFirst = FALSE, params = list(), mocks = function(rec) list()) {
  rec <- new.env(); rec$calls <- list()
  sim <- toySim(params)
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runDEoptim = function(...) {
      args <- list(...)
      rec$lower <- args$lower; rec$upper <- args$upper
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
