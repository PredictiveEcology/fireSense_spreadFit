## Shared by test-debugEvent.R and test-objectiveSettings.R.

recordingRSWD <- function(rec, value) {
  f <- function() {
    rec$args <- as.list(environment())
    value
  }
  formals(f) <- formals(runSpreadWithoutDEoptim)
  f
}

## the arguments the debug event gives runSpreadWithoutDEoptim(), and those the fit gives runDEoptim()
debugArgs <- function(params = list(), objects = list()) {
  rec <- new.env()
  sim <- toySim(c(list(mode = "debug"), params), objects)
  mockFitAndLedger(sim, new.env())
  mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, NULL))
  suppressMessages(SpaDES.core::spades(sim))
  rec$args
}
fitArgs <- function(params = list(), objects = list()) {
  rec <- new.env()
  sim <- toySim(params, objects)
  mockFitAndLedger(sim, rec)
  suppressMessages(SpaDES.core::spades(sim))
  rec$deArgs
}
