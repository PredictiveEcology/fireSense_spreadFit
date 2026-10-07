## After every fit, `postFitDiagnostics` turns what fireSenseUtils::runDEoptim() computed on the
## workers -- the re-score, the one-at-a-time profile, the uncapped simulations -- into sim objects,
## so the validation work of September 2026 is recovered from the simList, not from ad hoc scripts.
## mode "validate" adds `crossValidate`: two more fits, each predicting the years the other used.

m <- function(sim) SpaDES.core::params(sim)[[moduleName]]

## toyDE() with the attributes runDEoptim() attaches
toyDEWithDiagnostics <- function(nPar) {
  DE <- toyDE(nPar)
  attr(DE, "finalRescore") <- data.table::data.table(member = rep(1:7, each = 2), rep = 1:2,
                                                     value = rep(7:1, each = 2) + c(-0.5, 0.5))
  attr(DE, "profile") <- data.table::data.table(
    coef = c("", "CMDsm", "class1"), at = c(NA, 0, 0), isZero = c(FALSE, TRUE, TRUE),
    mean = c(1, 9, 1.1), delta = c(0, 8, 0.1), deltaSE = c(0, 1, 0.2))
  sims <- data.table::data.table(member = 1L, yr = "year2001", rep = 1:2, ids = 11L, size = 5, sim = c(5L, 10L))
  data.table::setattr(sims, "spreadProb",
                      list(fireSenseUtils::spreadProbSummary(c(0.1, 0.2, 0.25, 0.25), ceiling = 0.25)))
  attr(DE, "fitSims") <- sims
  DE
}

fittedWith <- function(params = list(), DEfun = toyDEWithDiagnostics) {
  rec <- new.env()
  sim <- toySim(c(list(stopIfNoPreRunFit = FALSE), params))
  mockFitAndLedger(sim, rec)
  mockInModule(sim, runDEoptim = function(...) {
    rec$deArgs <- list(...)
    DEfun(length(rec$deArgs$lower))
  })
  list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
}

test_that("the fit gets the chosen likelihood and asks for the diagnostics", {
  a <- fittedWith()$rec$deArgs
  expect_identical(a$sizeLik, "t")
  expect_identical(a$sizeLikDf, 5)
  expect_identical(a$weighted, FALSE)
  expect_identical(a$adWeight, "auto")
  expect_null(a$link)                                  # logistic3p: the objective's default
  expect_identical(a$profileReps, 10L)
  expect_identical(a$simulateMembers, 10L)
})

test_that("link 'logistic3pUpper' adds upperTail1 as the 2nd parameter and names the link", {
  a <- fittedWith(list(link = "logistic3pUpper", upperTailBounds = c(-0.5, 0.8)))$rec$deArgs
  expect_identical(a$link, "logistic3pUpper")
  expect_identical(names(a$lower)[1:2], c("maxAsymptote", "upperTail1"))
  expect_identical(names(a$upper), names(a$lower))
  expect_identical(unname(c(a$lower["upperTail1"], a$upper["upperTail1"])), c(-0.5, 0.8))
})

test_that("postFitDiagnostics makes the re-score, identifiability, sizes, saturation and convergence", {
  sim <- fittedWith()$sim
  ## re-score: member k's mean is 8 - k, one row per member, with the population's parameters
  expect_identical(nrow(sim$spreadFitRescore), 7L)
  expect_equal(sim$spreadFitRescore$reMean, 7:1)
  expect_equal(sim$spreadFitRescore$maxAsymptote, 1:7)
  ## identifiability: every covariate; the toy population holds them all at 0, so none is pinned,
  ## and none is identified even where dropping it matters (CMDsm)
  id <- sim$spreadFitIdentifiability
  ## yearSpreadSD too: "dropping" it sets the random effect's sd to 0, i.e. asks whether the effect matters
  expect_identical(id$coef, c("CMDsm", "youngAge", "class1", "class2", "nf", "yearSpreadSD"))
  expect_equal(id[coef == "CMDsm"]$dropMatters, TRUE)
  expect_equal(id[coef == "class1"]$dropMatters, FALSE)
  expect_false(any(id$identified))
  expect_identical(sim$spreadFitProfile$coef, c("", "CMDsm", "class1"))
  ## in-sample sizes and saturation, from the uncapped simulations
  expect_identical(nrow(sim$spreadFitSizes), 1L)
  expect_equal(sim$spreadFitSizes$totalAreaRatio, 1.5)
  expect_equal(sim$spreadFitLinkSaturation$propAtCeiling, 0.5)
  ## convergence: one row per block of the DEoptim result, in run order
  expect_equal(sim$spreadFitConvergence$bestval, c(30, 10, 20))
})

test_that("without the diagnostics' attributes, only what can be computed is made", {
  sim <- fittedWith(DEfun = toyDE)$sim
  expect_null(sim$spreadFitRescore)
  expect_null(sim$spreadFitProfile)
  expect_null(sim$spreadFitSizes)
  expect_false("identified" %in% names(sim$spreadFitIdentifiability))
  expect_identical(nrow(sim$spreadFitConvergence), 3L)
})

test_that("crossValidate fits each half of the years and predicts the other half, without the ledger", {
  rec <- new.env()
  rec$fits <- list(); rec$sims <- list()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE, mode = "validate", simulateMembers = 3L),
                list(studyAreaWithSpreadParams = toyLedger("9.9")))
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
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_identical(names(rec$fits), c("toyRun_cvFold1", "toyRun_cvFold2"))
  expect_identical(rec$fits$toyRun_cvFold1$years, "year2002")
  expect_identical(rec$fits$toyRun_cvFold2$years, "year2001")
  expect_identical(rec$fits$toyRun_cvFold1$profileReps, 0L)   # no profile of a validation fit
  ## each fold predicts the years its fit did not see, with the simulateMembers best members
  expect_identical(vapply(rec$sims, `[[`, "", "years"), c("year2001", "year2002"))
  expect_identical(nrow(rec$sims[[1]]$pop), 3L)
  ho <- sim$spreadFitHeldOut
  expect_identical(ho$sims$fold, 1:2)
  expect_equal(ho$score$fireBias, log10(2))
  expect_null(rec$geoArgs)                             # the ledger is not touched
  ## the simList is discarded in batch runs stopped after crossValidate, so the result must also
  ## reach disk (fireSense_spreadFit.R crossValidate event)
  heldOutPath <- file.path(toyFitDir(sim), "spreadFitHeldOut_toyRun.rds")
  expect_true(file.exists(heldOutPath))
  expect_equal(readRDS(heldOutPath), ho)
})

test_that("cvFolds() alternates years in year order", {
  expect_identical(cvFolds(c("year2003", "year2001", "year2002", "year2004")), c(1L, 1L, 2L, 2L))
})
