## Shared by test-climateOutside.R and test-powerLaw.R: a fit of the toy inputs under a spreadLink, and the Cache
## omitArgs the fit and the threshold calibration are given.

formCO <- "~ 0 + CMDsm + youngAge + class1 + class2 + nf"
refCol <- fireSenseUtils::spreadFitClimateRefTxt
linkCol <- fireSenseUtils::spreadLinkTxt

## the toy non-annual table named by data year, as the objective's covariate code needs (see test-spreadIntercept.R)
toyObjectsCO <- function(formula = formCO) {
  nonAnn <- toyObjects()$fireSense_nonAnnualSpreadFitCovariates
  list(fireSense_spreadFormula = formula,
       fireSense_nonAnnualSpreadFitCovariates = stats::setNames(nonAnn, "2000"))
}

## a fit on toy inputs; `rec` records what leaves the machine
fittedCO <- function(params = list(), formula = formCO, rec = new.env()) {
  sim <- toySim(params, toyObjectsCO(formula))
  mockFitAndLedger(sim, rec)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  list(sim = sim, rec = rec)
}

## what Cache(omitArgs =) the fit and the threshold calibration are given, and what the calibration receives
omitArgsCO <- function(params = list()) {
  rec <- new.env()
  sim <- toySim(c(list(SNLL_FS_thresh = NULL, iterThresh = 3L), params), toyObjectsCO())
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

