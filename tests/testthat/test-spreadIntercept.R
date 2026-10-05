## An optional free intercept (formula "~ 1 + ...", built by fireSense_dataPrepFit's `spreadIntercept`) with
## centred covariates. With an intercept: wide symmetric bounds for it, the covariate means computed once from
## the data the objective uses and handed to the fit, the threshold, the held-out simulation and the ledger row.
## Without one (the default, "~ 0 + ..."): nothing changes -- no intercept bound, no centre, no ledger column,
## and the cache keys of the fit, the threshold and the held-out simulation are what they were.

formI <- "~ 0 + CMDsm + youngAge + class1 + class2 + nf"
formIcpt <- sub("0 \\+", "1 +", formI, fixed = FALSE)
icptName <- fireSenseUtils::spreadInterceptTxt

## The toy non-annual table is named "year2001_year2002", which the real objective could not match to a year
## (the module names them by their first year, "2001"). Finding the centre runs the objective's covariate
## code for real, so the intercept tests name it as the module does.
toyObjectsI <- function(formula) {
  nonAnn <- toyObjects()$fireSense_nonAnnualSpreadFitCovariates
  list(fireSense_spreadFormula = formula,
       fireSense_nonAnnualSpreadFitCovariates = stats::setNames(nonAnn, "2001"))
}

## a fit on toy inputs, with the formula as given; `rec` records what leaves the machine
fittedI <- function(formula, params = list(), rec = new.env()) {
  sim <- toySim(c(list(stopIfNoPreRunFit = FALSE), params), toyObjectsI(formula))
  mockFitAndLedger(sim, rec)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  list(sim = sim, rec = rec)
}

test_that("estimateSpreadParams(): an intercept is bounded +/- upperAndLower, first among the coefficients", {
  annual <- list(year2001 = data.table::data.table(pixelID = 1L, CMDsm = 1, youngAge = 0))
  up <- estimateSpreadParams("~ 1 + CMDsm + youngAge + nf", annual, "upper", upperAndLower = 9)
  lo <- estimateSpreadParams("~ 1 + CMDsm + youngAge + nf", annual, "lower", upperAndLower = 9)
  expect_identical(up, c(maxAsymptote = 0.276, `(Intercept)` = 9, CMDsm = 9, youngAge = 0, nf = 9))
  expect_identical(lo, c(maxAsymptote = 0.25, `(Intercept)` = -9, CMDsm = 0, youngAge = -9, nf = -9))
  ## with upperTail1 it still comes right after the logistic parameters, and yearSpreadSD stays last
  up2 <- estimateSpreadParams("~ 1 + nf", annual, "upper", 9, upperTailBounds = c(-1, 1), yearSpreadSDBounds = c(0, 1))
  expect_identical(names(up2), c("maxAsymptote", "upperTail1", "(Intercept)", "nf", "yearSpreadSD"))
  ## it is the name fireSenseUtils defines
  expect_identical(names(up)[2], icptName)
  ## no intercept in the formula, none in the bounds
  expect_false(icptName %in% names(estimateSpreadParams("~ 0 + CMDsm + nf", annual, "upper", 9)))
})

test_that("spreadFitPrepare: with an intercept the default bounds include it; the centre is the data's mean", {
  p <- SpaDES.core::params(runEvents(runEvents(
    toySim(list(stopIfNoPreRunFit = FALSE), toyObjectsI(formIcpt)), "init"),
    "spreadFitPrepare"))[[moduleName]]
  expect_identical(names(p$upper)[1:2], c("maxAsymptote", icptName))
  expect_identical(unname(c(p$lower[[icptName]], p$upper[[icptName]])), c(-50, 50))
  expect_identical(names(p$lower), names(p$upper))
})

test_that("with an intercept the fit gets the centre, computed once from the data the objective uses", {
  out <- fittedI(formIcpt)
  a <- out$rec$deArgs
  centre <- out$sim$covCentre_spread
  expect_named(centre, c("CMDsm", "youngAge", "class1", "class2", "nf"))
  ## by hand, for CMDsm (rescaled with the range the fit stores, which the module may fix) and youngAge (as is)
  cmd <- c(10, 20, 30, 40, 25, 15, 12, 22, 32)
  rng <- out$sim$covMinMax_spread$CMDsm
  expect_equal(centre$CMDsm, mean((cmd - rng[1]) / (rng[2] - rng[1])))
  expect_equal(centre$youngAge, mean(c(0, 1, 0, 0, 0, 1, 0, 0, 1)))
  ## and it is exactly the one runDEoptim() receives, found from the data it receives
  expect_identical(a$covCentre, centre)
  expect_equal(centre, fireSenseUtils::spreadCovCentre(a$annualDTx1000, a$nonAnnualDTx1000, a$formulaToFit,
                                                       a$covMinMax, a$mutuallyExclusive))
  expect_identical(names(a$lower)[1:2], c("maxAsymptote", icptName))
})

test_that("the intercept's centre is in the ledger row, with the intercept among the parameter sets", {
  out <- fittedI(formIcpt)
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_true(fireSenseUtils::spreadFitCovCentreTxt %in% names(row))
  expect_identical(row[[fireSenseUtils::spreadFitCovCentreTxt]][[1]], out$sim$covCentre_spread)
  expect_identical(names(row$params[[1]])[1:4], c("maxAsymptote", "hillSlope1", "inflectionPoint1", icptName))
  expect_identical(out$sim$studyAreaWithSpreadParams, row)
})

test_that("off, there is no centre, no ledger column, and the fit's cache key omits it", {
  out <- fittedI(formI)
  expect_null(out$sim$covCentre_spread)
  expect_null(out$rec$deArgs$covCentre)
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_false(fireSenseUtils::spreadFitCovCentreTxt %in% names(row))
  expect_identical(setdiff(names(row), c("polygonID", "geometry")), fireSenseUtils::spreadFitAdditionalColNamesTxt)
  expect_false(icptName %in% names(row$params[[1]]))
  expect_false(icptName %in% names(out$rec$deArgs$lower))
})

## what Cache(omitArgs =) the fit and the threshold calibration are given
omitArgsOf <- function(formula, params = list()) {
  rec <- new.env()
  sim <- toySim(c(list(stopIfNoPreRunFit = FALSE), params), toyObjectsI(formula))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runSpreadWithoutDEoptim = function(...) { rec$rswd <- list(...); 777 },
    Cache = function(FUN, ..., omitArgs = NULL) {
      f <- if (is.function(FUN)) FUN(..., omitArgs = NULL) else FUN
      rec$omit <- c(rec$omit, list(omitArgs))
      f
    })
  suppressMessages(SpaDES.core::spades(sim))
  rec
}

test_that("the fit's and the threshold's cache keys omit covCentre only when it is NULL", {
  off <- omitArgsOf(formI, list(SNLL_FS_thresh = NULL, iterThresh = 3L))
  expect_true(all(vapply(off$omit, function(o) "covCentre" %in% o, logical(1))))
  expect_length(off$omit, 2L)                 # the threshold, then the fit
  on <- omitArgsOf(formIcpt, list(SNLL_FS_thresh = NULL, iterThresh = 3L))
  expect_false(any(vapply(on$omit, function(o) "covCentre" %in% o, logical(1))))
})

test_that("the threshold calibration is given the centre", {
  on <- omitArgsOf(formIcpt, list(SNLL_FS_thresh = NULL, iterThresh = 3L))
  expect_named(on$rswd$covCentre, c("CMDsm", "youngAge", "class1", "class2", "nf"))
  off <- omitArgsOf(formI, list(SNLL_FS_thresh = NULL, iterThresh = 3L))
  expect_null(off$rswd$covCentre)
})

test_that("runSpreadWithoutDEoptim() hands covCentre to the objective and to the gates", {
  src <- paste(deparse(body(runSpreadWithoutDEoptim)), collapse = "\n")
  ## the debug objective call, the gates, and the trials
  expect_equal(lengths(regmatches(src, gregexpr("covCentre = covCentre", src))), 3L)
  expect_true("covCentre" %in% names(formals(runSpreadWithoutDEoptim)))
  expect_null(formals(runSpreadWithoutDEoptim)$covCentre)
})

test_that("the held-out simulation and the fold's ledger row use the fit's centre", {
  run <- function(formula) {
    rec <- new.env(); rec$fnArgs <- list()
    sim <- toySim(list(heldOutFold = 1L, simulateMembers = 2L, stopIfNoPreRunFit = TRUE),
                  toyObjectsI(formula))
    mockFitAndLedger(sim, rec)
    mockInModule(sim, Cache = function(FUN, pop, fnArgs, ...) {
      if (!is.function(FUN)) return(FUN)
      rec$fnArgs <- fnArgs
      data.table::data.table(member = 1L, yr = "year2001", rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
    list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
  }
  on <- run(formIcpt)
  expect_identical(on$rec$fnArgs$covCentre, on$sim$covCentre_spread)
  expect_identical(on$sim$spreadFitHeldOut$fit[[fireSenseUtils::spreadFitCovCentreTxt]][[1]],
                   on$sim$covCentre_spread)
  off <- run(formI)
  expect_false("covCentre" %in% names(off$rec$fnArgs))
  expect_false(fireSenseUtils::spreadFitCovCentreTxt %in% names(off$sim$spreadFitHeldOut$fit))
})
