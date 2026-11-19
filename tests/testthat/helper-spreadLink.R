## Shared by test-climateOutside.R and test-powerLaw.R: a fit of the toy inputs under a spreadLink, and the Cache
## omitArgs the fit is given.

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

## What the module's Cache() does with the powerLaw calibration in a test that mocks Cache(): run it, without the
## Cache() arguments
runCalibration <- function(FUN, ...) {
  a <- list(...)
  do.call(FUN, a[setdiff(names(a), c("omitArgs", ".cacheExtra", "useCache"))])
}

## what Cache(omitArgs =) the fit is given
omitArgsCO <- function(params = list()) {
  rec <- new.env()
  sim <- toySim(params, toyObjectsCO())
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    Cache = function(FUN, ..., omitArgs = NULL) {
      if (identical(FUN, fireSenseUtils::calibratePowerLaw)) return(runCalibration(FUN, ..., omitArgs = omitArgs))
      f <- if (is.function(FUN)) FUN(..., omitArgs = NULL) else FUN
      rec$omit <- c(rec$omit, list(omitArgs))
      f
    })
  suppressMessages(SpaDES.core::spades(sim))
  rec
}

## the mean of the rescaled CMDsm over the 9 toy pixel-years, with the range the fit stores
cmdMean <- function(sim) {
  cmd <- c(10, 20, 30, 40, 25, 15, 12, 22, 32)
  rng <- sim$covMinMax_spread$CMDsm
  mean((cmd - rng[1]) / (rng[2] - rng[1]))
}

## records the validation data calls and the Plots() calls instead of simulating and drawing
mockValidation <- function(sim, rec) {
  rec$data <- list(); rec$plots <- list()
  mockInModule(sim,
    spreadFitValidationData = function(par, seed, ...) {
      a <- list(...)
      d <- data.table::data.table(year = names(a$historicalFires)[1], pixelID = 1L, observed = 1L,
                                  simulated = 0.5, p = 0.2)
      rec$data[[length(rec$data) + 1L]] <- list(par = par, seed = seed, args = a, d = d)
      d
    },
    ## spreadFitPrep() also calls Plots() (covariate histograms); only these figures are recorded
    Plots = function(data, fn, filename, ...) {
      if (!startsWith(filename, "spreadFit")) return(invisible(NULL))
      a <- list(...)
      rec$plots[[length(rec$plots) + 1L]] <- list(data = data, fn = fn, types = a$types, path = a$path,
                                                  filename = filename)
      invisible(NULL)
    })
}

fitWithPlots <- function(params = list(), objects = list()) {
  rec <- new.env()
  sim <- toySim(c(list(profileReps = 0L, simulateMembers = 0L), params), objects)
  mockFitAndLedger(sim, rec)
  mockValidation(sim, rec)
  list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
}
