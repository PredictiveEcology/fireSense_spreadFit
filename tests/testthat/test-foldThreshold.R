## A held-out fold's fit used the SNLL threshold calibrated on ALL years. The threshold bounds the SNLL
## of the two largest fire years of the data it is used on, and a fold's two largest years are not the
## full data's: on 2026-09-29 the fits for ELF 4.3 fold 2 and 5.2.1 fold 1 never passed it, returned the
## fail value 1e6 for all 5000 generations, and still wrote a held-out score from random parameters.

foldRun <- function(k, popval = NULL, runName = "toyRun") {
  rec <- new.env()
  rec$calibratedOn <- list(); rec$fitThresh <- list()
  sim <- toySim(list(heldOutFold = k, SNLL_FS_thresh = NULL, simulateMembers = 3L),
                objects = list(.runName = runName))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    ## the calibrated threshold records, and depends on, the years it was calibrated on
    runSpreadWithoutDEoptim = function(historicalFires, ...) {
      rec$calibratedOn[[length(rec$calibratedOn) + 1L]] <- names(historicalFires)
      100 * length(historicalFires)
    },
    runDEoptim = function(...) {
      a <- list(...)
      rec$fitThresh[[a$runName]] <- as.vector(a$thresh)
      de <- toyDE(length(a$lower))
      if (!is.null(popval)) de[[length(de)]]$member$popval[] <- popval
      de
    },
    Cache = function(FUN, pop, fnArgs, ...) {
      if (!is.function(FUN)) return(FUN)
      yr <- names(fnArgs$historicalFires)
      data.table::data.table(member = 1L, yr = yr, rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
  list(sim = sim, rec = rec)
}

test_that("a held-out fold's fit uses a threshold calibrated on the years it fits", {
  out <- foldRun(1L)
  suppressMessages(SpaDES.core::spades(out$sim))
  rec <- out$rec
  ## heldOutFold 1 fits year2002 only (see test-heldOutFold.R)
  expect_true(list("year2002") %in% rec$calibratedOn)
  expect_identical(rec$fitThresh$toyRun_cvFold1, 100)
})

test_that("a fit whose final population is all fail values stops instead of scoring random parameters", {
  out <- foldRun(1L, popval = 1e6, runName = "toyFailedFit")
  expect_error(suppressMessages(SpaDES.core::spades(out$sim)), "fail value")
  heldOutPath <- file.path(SpaDES.core::outputPath(out$sim), moduleName, "spreadFitHeldOut_toyFailedFit_fold1.rds")
  expect_false(file.exists(heldOutPath))
})
