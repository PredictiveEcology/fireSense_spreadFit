## A project passes the run name as a PARAMETER (`.runName` in `params$.globals`), the way
## Biomass_core receives it. The module used to read only the input object, so every
## production fit had a NULL name: progress files were `DEoptimProgress_.csv` and clusters
## recorded worker memory under no name.

test_that("the run name comes from the parameter when there is no object", {
  sim <- toySim(list(.runName = "x"), list(.runName = NULL))
  expect_null(sim$.runName)
  expect_identical(runNameOf(sim), "x")
})

test_that("the object is used when the parameter is not set", {
  expect_identical(runNameOf(toySim()), "toyRun")
})

test_that("the parameter wins over the object", {
  expect_identical(runNameOf(toySim(list(.runName = "x"))), "x")
})

test_that("a fit receives the parameter's name, and names its progress file with it", {
  rec <- new.env()
  sim <- toySim(list(simulateMembers = 3L, stopIfNoPreRunFit = FALSE, .runName = "fromParam"), list(.runName = NULL))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runDEoptim = function(...) {
      a <- list(...)
      rec$runName <- a$runName
      rec$progressFile <- a$progressFile
      toyDE(length(a$lower))
    },
    Cache = function(FUN, ..., omitArgs = NULL) {
      if (!is.function(FUN)) return(FUN)
      data.table::data.table(member = 1L, yr = "year2001", rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
  out <- suppressMessages(SpaDES.core::spades(sim))
  expect_identical(rec$runName, "fromParam")
  expect_identical(basename(rec$progressFile), "DEoptimProgress_fromParam.csv")
})
