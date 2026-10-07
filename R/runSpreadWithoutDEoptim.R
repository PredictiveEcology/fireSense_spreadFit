#' Evaluate the spread objective function without DEoptim
#'
#' In "debug" mode, runs `fireSenseUtils::.objfunSpreadFit()` once per parameter set, with plots.
#' Otherwise calibrates the SNLL threshold: evaluates `iterThresh` random parameter sets with no early
#' stop, in forked processes (drawn from the logistic's active range, see `drawActivePars()`), and returns `pickThreshold()` of their first-block SNLLs: `thresholdMargin`
#' times the best trial that did not saturate.
#'
#' @param iterThresh integer; number of random parameter sets.
#' @param thresholdMargin numeric; the threshold is this multiple of the best usable trial's first-block
#'   average annual SNLL (the unit of `thresh` in the objective).
#' @param lower,upper numeric; bounds the random parameter sets are drawn between.
#' @param maxDrawsPerTrial integer; the calibration keeps only draws that pass the objective's spreadProb
#'   gates in the first block of years (`fireSenseUtils::spreadProbGates()`), until `iterThresh` are
#'   kept or `maxDrawsPerTrial * iterThresh` have been drawn. A constant of the function, not a module
#'   parameter. With none kept the result is `NA` (no early stop, see `noEarlyStopThreshold()`).
#' @param fireSense_spreadFormula character; passed to `FS_formula`.
#' @param flammableRTM `SpatRaster`; passed to `landscape`.
#' @param annualDTx1000,nonAnnualDTx1000,fireBufferedListDT,historicalFires elements of the list made
#'   by `fireSenseUtils::covsX1000AndSetDF()`.
#' @param mutuallyExclusive named list of mutually exclusive covariates.
#' @param doObjFunAssertions logical; passed to `doAssertions`.
#' @param covMinMax `data.table` of covariate min and max, or NULL for no rescaling.
#' @param covCentre named list of covariate centres (`fireSenseUtils::spreadCovCentre()`) subtracted from the
#'   rescaled covariates, as in the fit; `NULL` (no intercept) centres nothing. Passed to `.objfunSpreadFit()`
#'   and `fireSenseUtils::spreadProbGates()`, so the threshold is calibrated on the fit's objective.
#' @param objfunFireReps integer; passed to `Nreps`, the replicates per fire.
#' @param maxFireSpread numeric; upper limit on mean spread probability.
#' @param weighted logical; weight the SNLL by log fire size. Also used in the rough threshold estimate.
#' @param tests character; objective function tests, e.g. "SNLL_FS", "adTest".
#' @param formulaToFit character; the spread formula.
#' @param pars optional numeric vector, or list of them, to evaluate instead of random sets.
#' @param plot.it passed to `.objfunSpreadFit()` in "debug" mode.
#' @param mode character; if it includes "debug", the debug branch runs.
#' @param seed integer or NULL; `set.seed()` value. NULL draws one at random.
#' @param escapeSizeHa passed to `.objfunSpreadFit()`: fit escaped fires only (NULL: any fire over 1 pixel).
#' @param sizeLik,sizeLikDf,adWeight,link,jumpTries,jumpMeanDist,yearAreaWeight,areaDistWeight,penaliseRunaways,runawayEdgeFrac,runawayEdgeMin passed to
#'   `.objfunSpreadFit()`, so the threshold is calibrated on the objective the fit uses. The defaults are
#'   that function's.
#' @return the calibrated threshold (numeric, or NA if no trial was usable); NULL in "debug" mode.
runSpreadWithoutDEoptim <- function(iterThresh, lower, upper, fireSense_spreadFormula, flammableRTM,
                                    annualDTx1000, nonAnnualDTx1000, fireBufferedListDT,
                                    mutuallyExclusive = list("youngAge" = "vegPC"),
                                    doObjFunAssertions = getOption("fireSenseUtils.assertions", TRUE),
                                    historicalFires, covMinMax, covCentre = NULL, objfunFireReps, maxFireSpread,
                                    weighted = TRUE, tests = c("snll_fs", "adtest"),
                                    formulaToFit,
                                    pars = NULL, plot.it = TRUE, mode = "fit",
                                    seed = NULL, escapeSizeHa = NULL,
                                    sizeLik = "kde", sizeLikDf = 5, adWeight = "auto", link = NULL,
                                    jumpTries = 0, jumpMeanDist = 0, yearAreaWeight = 0, areaDistWeight = 0,
                                    penaliseRunaways = TRUE,
                                    runawayEdgeFrac = fireSenseUtils::fireSenseRunawayEdgeFrac,
                                    runawayEdgeMin = fireSenseUtils::fireSenseRunawayEdgeMin,
                                    thresholdMargin = 2, maxDrawsPerTrial = 100) {
  ## The threshold this returns becomes `thresh` in runDEoptim(), so it is part of every cached
  ## DEoptim generation's key. With a seed drawn here, a single cache miss on the estimateThreshold
  ## event re-drew the threshold and invalidated EVERY cached generation for that ELF: on 2026-09-16
  ## a restarted fit went 1236 -> 1416 and recomputed from generation 1, losing ~17 h, while one that
  ## hit the cache kept 1705 and replayed 797 generations in ~35 min. A caller that passes a seed
  ## derived from the ELF gets a reproducible threshold, so a miss costs only this estimate.
  ## NULL keeps the old behaviour for callers that do not care (e.g. the module's `debug` event).
  if (is.null(seed)) seed <- sample(1e6, 1)
  set.seed(seed)

  n <- iterThresh ## the more you do, the lower the resulting threshold

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
  drawFromActiveRange <- is.null(pars) && !"debug" %in% mode
  if (is.null(pars)) {
    ## do NOT re-draw here: that discarded the seed set above, which is what made the threshold
    ## irreproducible even when the caller asked for a specific seed
    print(paste("seed used for runSpreadWithoutDEoptim is ", seed))
    ## the fit branch draws below, from the logistic's active range; "debug" evaluates uniform draws
    if (!drawFromActiveRange) pars <- lapply(1:n, function(x) runif(length(lower), lower, upper))
  }
  if (!drawFromActiveRange) pars <- namePars(pars)


  if ("debug" %in% mode) {
    a <- list()
    for (i in seq(pars)) {
      print(paste(i, "logit params:", paste(round(pars[[i]], 2), collapse = ", ")))
      a[[i]] <- .objfunSpreadFit(par = pars[[i]],
                                 thresh = Inf,
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
                                 covMinMax = covMinMax, covCentre = covCentre,
                                 Nreps = objfunFireReps,
                                 maxFireSpread = maxFireSpread,
                                 verbose = TRUE,
                                 weighted = weighted,
                                 escapeSizeHa = escapeSizeHa,
                                 sizeLik = sizeLik, sizeLikDf = sizeLikDf, adWeight = adWeight, link = link,
                                 jumpTries = jumpTries, jumpMeanDist = jumpMeanDist,
                                 yearAreaWeight = yearAreaWeight, areaDistWeight = areaDistWeight,
                                 penaliseRunaways = penaliseRunaways,
                                 runawayEdgeFrac = runawayEdgeFrac, runawayEdgeMin = runawayEdgeMin,
                                 plot.it = plot.it
      )
    }
  } else {
    message("SNLL_FS_thresh not specified. Self calibrating threshold value for runDEoptim (n=", n, ")")

    # Check for being in a future
    a <- future::plan()
    if (is(a, "FutureStrategy") && !is(a, "sequential")) {
      coresToUse <- nbrOfWorkers()
    } else {
      activeThreads <- clusters::numActiveThreads()
      detCores <- detectCores()
      # Each fork is a copy-on-write image of this process that diverges as the
      # garbage collector marks the heap, so budget one heap per fork; the gc()
      # call also leaves less garbage to copy. See thresholdForks().
      heapMB <- sum(gc()[, 2])
      availMB <- pemisc::availableMemory()
      availMB <- if (length(availMB)) availMB / 1e6 else NULL
      coresToUse <- thresholdForks(heapMB = heapMB, availMB = availMB,
                                   nPars = if (drawFromActiveRange) n else length(pars),
                                   detCores = detCores, activeThreads = activeThreads,
                                   mcCores = getOption("mc.cores"))
      withr::local_options("mc.cores" = coresToUse)
    }
    message("Using ", coresToUse, " cores",
            if (exists("heapMB", inherits = FALSE))
              paste0(" (R heap ", round(heapMB / 1024), " GB; host has ",
                     if (length(availMB)) paste(round(availMB / 1024), "GB available") else "unknown memory", ")"),
            ".")

    if (drawFromActiveRange) {
      maxDraws <- as.integer(maxDrawsPerTrial * n)
      gates <- function(parList) {
        spreadProbGates(parList, annualDTx1000 = annualDTx1000, nonAnnualDTx1000 = nonAnnualDTx1000,
                        historicalFires = historicalFires, formulaToFit = formulaToFit,
                        covMinMax = covMinMax, covCentre = covCentre, mutuallyExclusive = mutuallyExclusive,
                        maxFireSpread = maxFireSpread, escapeSizeHa = escapeSizeHa,
                        landscape = flammableRTM, link = link, doAssertions = FALSE)$pass
      }
      active <- drawActivePars(n, lower, upper, maxDraws, gates, cores = coresToUse)
      message("  accepted ", active$accepted, " of ", active$drawn, " draws (cap ", maxDraws, " = ",
              maxDrawsPerTrial, " x iterThresh) that pass the spreadProb gates of the first block of years")
      if (!active$accepted) {
        message("  no draw passed the spreadProb gates; no threshold calibrated")
        return(NA_real_)
      }
      pars <- namePars(active$pars)
    }

    st1 <- system.time({
      objSpreadFit <- mcmapply(mc.cores = coresToUse,
                    mc.preschedule = FALSE, SIMPLIFY = FALSE,
                    par = pars, FUN = trialFirstBlock,
                    MoreArgs = list(
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
                      covMinMax = covMinMax, covCentre = covCentre,
                      Nreps = objfunFireReps,
                      maxFireSpread = maxFireSpread,
                      weighted = weighted, escapeSizeHa = escapeSizeHa,
                      sizeLik = sizeLik, sizeLikDf = sizeLikDf, adWeight = adWeight, link = link,
                      jumpTries = jumpTries, jumpMeanDist = jumpMeanDist,
                      yearAreaWeight = yearAreaWeight, areaDistWeight = areaDistWeight,
                      penaliseRunaways = penaliseRunaways,
                      plot.it = FALSE)
      )
    })

    first <- lapply(objSpreadFit, function(x) if (is.list(x)) x else list(annual = NA_real_, saturated = NA))
    annual <- vapply(first, function(x) as.numeric(x$annual)[1], numeric(1))
    saturated <- vapply(first, function(x) isTRUE(x$saturated), logical(1))
    threshToUse <- suppressWarnings(pickThreshold(annual, saturated, margin = thresholdMargin))
    nUsable <- sum(is.finite(annual) & !saturated)
    message("  ", nUsable, " of ", length(pars), " trials usable (first-block avg annual SNLL, no early stop: ",
            if (nUsable) paste0("best ", min(annual[is.finite(annual) & !saturated])) else "none",
            "); using SNLL_FS_thresh value: ", threshToUse, " (thresholdMargin ", thresholdMargin, ")")
    return(threshToUse)
  }
}
