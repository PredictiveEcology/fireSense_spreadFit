#' Evaluate the spread objective function without DEoptim
#'
#' The module's `debug` event: runs `fireSenseUtils::.objfunSpreadFit()` once per parameter set, with plots.
#'
#' @param nSets integer; number of random parameter sets, drawn uniformly between `lower` and `upper`.
#'   Not used with `pars`.
#' @param lower,upper numeric; bounds the random parameter sets are drawn between.
#' @param fireSense_spreadFormula character; passed to `FS_formula`.
#' @param flammableRTM `SpatRaster`; passed to `landscape`.
#' @param annualDTx1000,nonAnnualDTx1000,fireBufferedListDT,historicalFires elements of the list made
#'   by `fireSenseUtils::covsX1000AndSetDF()`.
#' @param mutuallyExclusive named list of mutually exclusive covariates.
#' @param doObjFunAssertions logical, no default: the module passes its parameter of that name, so one
#'   place (`options(fireSenseUtils.assertions = )`) sets it. Passed to `doAssertions`.
#' @param covMinMax `data.table` of covariate min and max, or NULL for no rescaling.
#' @param covCentre named list of covariate centres (`fireSenseUtils::spreadCovCentre()`) subtracted from the
#'   rescaled covariates, as in the fit; `NULL` (no intercept) centres nothing. Passed to `.objfunSpreadFit()`.
#' @param climateRef named list, the fit's climate reference (`fireSenseUtils::spreadClimateRef()`; `spreadLink =
#'   "climateOutside"`), or `NULL` (the default link). Passed to `.objfunSpreadFit()` like `covCentre`.
#' @param powerLawGamma numeric; the `"powerLaw"` link's `gamma` (`fireSenseUtils::powerLawLink()`), used with
#'   `link = "powerLaw"`. Passed to `.objfunSpreadFit()`.
#' @param powerLawShape numeric; the `"powerLaw"` link's `shape`, passed like `powerLawGamma`.
#' @param lowerAsymptote,richardsNu numeric; the lower asymptote of the `"powerLawSmooth"` and `"richards"` links and the
#'   `nu` of `"richards"` (`fireSenseUtils::powerLawSmoothLink()`, `richardsLink()`), passed like `powerLawGamma`.
#' @param objfunFireReps integer; passed to `Nreps`, the replicates per fire.
#' @param maxFireSpread numeric; upper limit on mean spread probability.
#' @param weighted passed to `.objfunSpreadFit()`.
#' @param tests character; objective function tests, e.g. "crps", "adTest".
#' @param formulaToFit character; the spread formula.
#' @param pars optional numeric vector, or list of them, to evaluate instead of random sets.
#' @param plot.it passed to `.objfunSpreadFit()`.
#' @param seed integer or NULL; `set.seed()` value. NULL draws one at random.
#' @param escapeSizeHa passed to `.objfunSpreadFit()`: fit escaped fires only (NULL: any fire over 1 pixel).
#' @param sizeLik,sizeLikDf,adWeight,link,jumpTries,jumpMeanDist,yearAreaWeight,areaDistWeight,penaliseRunaways,runawayEdgeFrac,runawayEdgeMin,runawayBufferMultiple passed to
#'   `.objfunSpreadFit()`, so it runs the objective the fit uses. The defaults are that function's.
#' @param wholeLandscape passed to `.objfunSpreadFit()`: spread every fire on all pixels of the tables.
#' @return `NULL`.
runSpreadWithoutDEoptim <- function(nSets = 1L, lower, upper, fireSense_spreadFormula, flammableRTM,
                                    annualDTx1000, nonAnnualDTx1000, fireBufferedListDT,
                                    mutuallyExclusive = list("youngAge" = "vegPC"),
                                    doObjFunAssertions,
                                    historicalFires, covMinMax, covCentre = NULL, climateRef = NULL, objfunFireReps, maxFireSpread,
                                    weighted = TRUE, tests = "crps",
                                    formulaToFit,
                                    pars = NULL, plot.it = TRUE,
                                    seed = NULL, escapeSizeHa = NULL,
                                    sizeLik = "kde", sizeLikDf = 5, adWeight = "auto", link = NULL,
                                    powerLawGamma = fireSenseUtils::spreadPowerLawGamma,
                                    powerLawShape = fireSenseUtils::spreadPowerLawShape,
                                    lowerAsymptote = fireSenseUtils::spreadProbFloor,
                                    richardsNu = fireSenseUtils::spreadRichardsNu, jumpTries = 0, jumpMeanDist = 0, yearAreaWeight = 0, areaDistWeight = 0,
                                    crpsWeights = NULL, gateRefusal = "minLik",
                                    penaliseRunaways = TRUE,
                                    runawayEdgeFrac = fireSenseUtils::fireSenseRunawayEdgeFrac,
                                    runawayEdgeMin = fireSenseUtils::fireSenseRunawayEdgeMin,
                                    runawayBufferMultiple = NULL, wholeLandscape = FALSE) {
  if (is.null(seed)) seed <- sample(1e6, 1)
  set.seed(seed)

  ## .objfunSpreadFit (fireSenseUtils) tells yearSpreadSD apart from a logistic parameter only by
  ## name (it is always the trailing element of lower/upper); name unnamed sets that match `lower`'s
  ## length, whether drawn above or supplied via `pars`, or every trial with yearSpreadSD errors.
  namePars <- function(pars) {
    if (!is(pars, "list")) pars <- list(pars)
    lapply(pars, function(p) {
      if (is.null(names(p)) && length(p) == length(lower)) names(p) <- names(lower)
      p
    })
  }
  if (is.null(pars)) {
    print(paste("seed used for runSpreadWithoutDEoptim is ", seed))
    pars <- lapply(seq_len(nSets), function(x) runif(length(lower), lower, upper))
  }
  pars <- namePars(pars)

  for (i in seq(pars)) {
    print(paste(i, "logit params:", paste(round(pars[[i]], 2), collapse = ", ")))
    .objfunSpreadFit(par = pars[[i]],
                     FS_formula = fireSense_spreadFormula,
                     landscape = flammableRTM,
                     annualDTx1000 = annualDTx1000,
                     nonAnnualDTx1000 = nonAnnualDTx1000,
                     fireBufferedListDT = fireBufferedListDT,
                     mutuallyExclusive = mutuallyExclusive,
                     doAssertions = doObjFunAssertions,
                     historicalFires = historicalFires,
                     formulaToFit = formulaToFit,
                     tests = tests,
                     covMinMax = covMinMax, covCentre = covCentre, climateRef = climateRef,
                     Nreps = objfunFireReps,
                     maxFireSpread = maxFireSpread,
                     verbose = TRUE,
                     weighted = weighted,
                     escapeSizeHa = escapeSizeHa,
                     sizeLik = sizeLik, sizeLikDf = sizeLikDf, adWeight = adWeight, link = link,
                     powerLawGamma = powerLawGamma, powerLawShape = powerLawShape,
                     lowerAsymptote = lowerAsymptote, richardsNu = richardsNu, jumpTries = jumpTries, jumpMeanDist = jumpMeanDist,
                     yearAreaWeight = yearAreaWeight, areaDistWeight = areaDistWeight,
                     crpsWeights = crpsWeights, gateRefusal = gateRefusal,
                     penaliseRunaways = penaliseRunaways,
                     runawayEdgeFrac = runawayEdgeFrac, runawayEdgeMin = runawayEdgeMin,
                     runawayBufferMultiple = runawayBufferMultiple, wholeLandscape = wholeLandscape,
                     plot.it = plot.it)
  }
  invisible(NULL)
}
