## A fit's per-generation progress goes to a CSV in the run's own output folder, one file per fit
## (the main fit and each cross-validation fold), without any setting. runDEoptim() hands it to
## clusters::DEoptimIterative().

progressFileRun <- function(params = list(), preRun = !is.null(params$heldOutFold)) {
  rec <- new.env()
  sim <- toySim(c(list(simulateMembers = 3L, stopIfNoPreRunFit = preRun), params))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runDEoptim = function(...) {
      a <- list(...)
      rec$progressFiles[[a$runName]] <- a$progressFile
      toyDE(length(a$lower))
    },
    Cache = function(FUN, ..., omitArgs = NULL) {
      if (!is.function(FUN)) { rec$omitArgs <- omitArgs; return(FUN) }
      data.table::data.table(member = 1L, yr = "year2001", rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
  rec$progressFiles <- list()
  list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
}

test_that("each fit gets its own progress file in the module's output folder", {
  out <- progressFileRun(list(heldOutFold = 1L))
  dir <- file.path(SpaDES.core::outputPath(out$sim), moduleName)
  expect_identical(out$rec$progressFiles$toyRun_cvFold1, file.path(dir, "DEoptimProgress_toyRun_cvFold1.csv"))
  out2 <- progressFileRun(list(heldOutFold = 2L))
  expect_identical(out2$rec$progressFiles$toyRun_cvFold2,
                   file.path(SpaDES.core::outputPath(out2$sim), moduleName, "DEoptimProgress_toyRun_cvFold2.csv"))
  expect_true(dir.exists(dir))
})

test_that("the main fit's progress file is named for the run, and is not in the cache key", {
  out <- progressFileRun()
  expect_identical(out$rec$progressFiles$toyRun,
                   file.path(SpaDES.core::outputPath(out$sim), moduleName, "DEoptimProgress_toyRun.csv"))
  expect_true("progressFile" %in% out$rec$omitArgs)
})
