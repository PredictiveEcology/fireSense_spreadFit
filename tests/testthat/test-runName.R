## A project passes the run name as a PARAMETER (`.runName` in `params$.globals`), the way
## Biomass_core receives it. The module used to read only the input object, so every
## production fit had a NULL name: progress files were `DEoptimProgress_.csv` and clusters
## recorded worker memory under no name.

test_that("an input object .runName is no longer read: the name is the parameter, NA when unset", {
  sim <- toySim(list(.runName = NULL), list(.runName = "fromObject"))
  expect_identical(SpaDES.core::P(sim, moduleName)$.runName, NA_character_)
})

test_that("a fit receives the parameter's name, and names its progress file with it", {
  rec <- new.env()
  sim <- toySim(list(simulateMembers = 3L, .runName = "fromParam"))
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
