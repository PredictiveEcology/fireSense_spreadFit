## `heldOutFold` runs ONE cross-validation fold as its own job, instead of both folds in the same
## job (mode "validate"'s `crossValidate`). See fireSense_spreadFit.R init and R/fitSpread.R
## crossValidateSpreadOneFold().

heldOutFoldRun <- function(k, params = list()) {
  rec <- new.env()
  rec$fits <- list(); rec$sims <- list()
  sim <- toySim(c(list(heldOutFold = k, simulateMembers = 3L, stopIfNoPreRunFit = TRUE), params))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runDEoptim = function(...) {
      a <- list(...)
      rec$fits[[a$runName]] <- list(years = names(a$historicalFires), profileReps = a$profileReps)
      toyDE(length(a$lower))
    },
    ## Cache() wraps both the fit, already evaluated here, and the held-out simulation
    Cache = function(FUN, pop, fnArgs, ...) {
      if (!is.function(FUN)) return(FUN)
      rec$sims[[length(rec$sims) + 1L]] <- list(pop = pop, years = names(fnArgs$historicalFires))
      yr <- names(fnArgs$historicalFires)
      data.table::data.table(member = 1L, yr = yr, rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
  list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
}

test_that("heldOutFold = 1 fits fold 2's years, scores only fold 1, writes a fold-specific file", {
  out <- heldOutFoldRun(1L)
  sim <- out$sim; rec <- out$rec
  ## exactly one fit is run, on the OTHER fold's years (fold 1 is year2001; see cvFolds() test)
  expect_identical(names(rec$fits), "toyRun_cvFold1")
  expect_identical(rec$fits$toyRun_cvFold1$years, "year2002")
  expect_identical(rec$fits$toyRun_cvFold1$profileReps, 0L)
  ho <- sim$spreadFitHeldOut
  expect_identical(ho$sims$fold, 1L)
  expect_identical(vapply(rec$sims, `[[`, "", "years"), "year2001")   # predicts fold 1's held-out year
  expect_null(rec$geoArgs)                                           # the ledger is never touched
  heldOutPath <- file.path(SpaDES.core::outputPath(sim), moduleName, "spreadFitHeldOut_toyRun_fold1.rds")
  expect_true(file.exists(heldOutPath))
  expect_equal(readRDS(heldOutPath), ho)
  ## no full-fit file, and no mode-"validate" (both-fold) file
  expect_false(file.exists(file.path(SpaDES.core::outputPath(sim), moduleName, "spreadFitHeldOut_toyRun.rds")))
})

test_that("heldOutFold = 2 fits fold 1's years, scores only fold 2, writes a fold-specific file", {
  out <- heldOutFoldRun(2L)
  sim <- out$sim; rec <- out$rec
  expect_identical(names(rec$fits), "toyRun_cvFold2")
  expect_identical(rec$fits$toyRun_cvFold2$years, "year2001")
  ho <- sim$spreadFitHeldOut
  expect_identical(ho$sims$fold, 2L)
  expect_identical(vapply(rec$sims, `[[`, "", "years"), "year2002")
  heldOutPath <- file.path(SpaDES.core::outputPath(sim), moduleName, "spreadFitHeldOut_toyRun_fold2.rds")
  expect_true(file.exists(heldOutPath))
  expect_equal(readRDS(heldOutPath), ho)
})

test_that("heldOutFold never schedules or runs the full fit ('run') or writes the ledger", {
  out <- heldOutFoldRun(1L)
  expect_false("run" %in% SpaDES.core::events(out$sim)$eventType)
  expect_null(out$rec$geoArgs)
})

test_that("the held-out object carries the fold's fit in the ledger row's structure", {
  out <- heldOutFoldRun(1L)
  fold <- out$sim$spreadFitHeldOut
  runSim <- toySim(list(stopIfNoPreRunFit = FALSE))      # the `run` event's ledger row
  mockFitAndLedger(runSim)
  run <- suppressMessages(SpaDES.core::spades(runSim))$studyAreaWithSpreadParams
  fit <- fold$fit
  expect_identical(names(fit), names(run))                # same columns, same order
  expect_identical(vapply(fit, function(x) class(x)[1], ""), vapply(run, function(x) class(x)[1], ""))
  expect_identical(lapply(fit$params, class), lapply(run$params, class))
  expect_identical(class(fit), class(run))
  expect_identical(sf::st_crs(fit), sf::st_crs(run))
  expect_identical(fit$polygonID, run$polygonID)
  ## all simulateMembers (3) members, parameter columns as the ledger has them
  expect_identical(names(fit$params[[1]]), names(run$params[[1]]))
  expect_identical(NROW(fit$params[[1]]), 3L)
  expect_length(fit$objFunVal[[1]], 3L)
  expect_identical(fit$covMinMax_spread[[1]], run$covMinMax_spread[[1]])
  expect_identical(fold$formula, "~ 0 + CMDsm + youngAge + class1 + class2 + nf")
  expect_identical(fold$link, "logistic3p")
  heldOutPath <- file.path(SpaDES.core::outputPath(out$sim), moduleName, "spreadFitHeldOut_toyRun_fold1.rds")
  expect_identical(names(readRDS(heldOutPath)), c("sims", "score", "fit", "heldOutFold", "fitYears", "heldOutYears",
                                                   "formula", "link"))
  expect_identical(fold$heldOutFold, 1L)
  expect_identical(fold$fitYears, "year2002")
  expect_identical(fold$heldOutYears, "year2001")
})
