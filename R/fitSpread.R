#' Fit the spread model
#'
#' The module's one `fireSenseUtils::runDEoptim()` call, shared by the `run` event (every year) and
#' `crossValidate` (one fold's years). After the fit, `runDEoptim()` re-scores the final population
#' and, with `profileReps`/`simulateMembers` above 0, profiles the best member and simulates the best
#' members' fires, on the same workers; `makeFitDiagnostics()` turns those into `sim` objects.
#'
#' @param sim a `simList`.
#' @param covs the covariates as integers x 1000 (`mod$covsX1000`, or a subset of its years).
#' @param thresh the objective's early-stop threshold, `mod$thresh`.
#' @param runName character; labels the run and its cache entry. Also names the fit's progress file,
#'   `DEoptimProgress_<runName>.csv` in the module's output folder.
#' @param diagnostics logical; `FALSE` skips the profile and the simulations.
#' @param passingPars matrix or `NULL`; the threshold calibration's parameter sets that pass `thresh`
#'   (`mod$thresholdPassingPars`). They are the first members of DEoptim's starting population, on the first
#'   attempt and on the retry (`withPassingFirst()`), unless the `initialpop` parameter is set or the emulator
#'   chose the start (`initialPopMethod = "emulator"`). The folds
#'   pass none: their first block is other years.
#' @return the `runDEoptim()` result. Stops if every member of the final population has the
#'   objective's fail value: no parameter set ever passed `thresh`, so there is no fit.
fitSpread <- function(sim, covs, thresh, runName, diagnostics = TRUE, passingPars = NULL) {
  stopifnot("`thresh` is NA: it must be a number, or Inf for no early stop" = !anyNA(thresh))
  if (!is.null(P(sim)$cores) && !any(is.na(P(sim)$cores)) &&
      identical(sort(unique(P(sim)$cores)), sort(P(sim)$cores))) {
    best <- list(cluster = P(sim)$cores)
  } else {
    best <- list(cluster = P(sim)$cores,
                 bestCluster = as.data.table(table(P(sim)$cores)))
  }
  messageDF(best$bestCluster)
  fnName <- paste0("runDEoptim_", runName, "_", P(sim)$.rep)
  runawayMultiple <- runawayBufferMultipleOrNULL(P(sim)$runawayBufferMultiple)
  ## one progress file per fit (the main fit, each fold), in the run's own output folder
  progressFile <- file.path(checkPath(file.path(outputPath(sim), currentModule(sim)), create = TRUE),
                            paste0("DEoptimProgress_", runName, ".csv"))
  ## runDEoptim()'s arguments, once: the fit gets them, and the emulator (`initialPopMethod = "emulator"`)
  ## evaluates the objective they define
  deArgs <- list(landscape = sim$rasterToMatch,
                   annualDTx1000 = covs$annualDTx1000,
                   nonAnnualDTx1000 = covs$nonAnnualDTx1000,
                   fireBufferedListDT = covs$fireBufferedListDT,
                   historicalFires = covs$historicalFires,
                   itermax = P(sim)$iterDEoptim,
                   iterStep = 1L, ## hard-coded: iterStep > 1 crashed fits (clusters' c = 0 adaptation, see NEWS)
                   ## the cluster's size is the population size; see the parameter's doc
                   nCoresNeeded = P(sim)$nCoresNeeded,
                   trace = P(sim)$trace,
                   initialpop = P(sim)$initialpop,
                   strategy = P(sim)$strategy,
                   cores = best$cluster,
                   doObjFunAssertions = P(sim)$doObjFunAssertions,
                   paths = getPaths(),
                   libPath = normPath(P(sim)$libPathDEoptim),
                   logPath = logPath(sim), ## TODO (#6): use tempdir()
                   lower = P(sim)$lower,
                   upper = P(sim)$upper,
                   mutuallyExclusive = P(sim)$mutuallyExclusiveCols, ## TODO: test
                   formulaToFit = sim$fireSense_spreadFormula,
                   covMinMax = sim$covMinMax_spread,
                   covCentre = sim$covCentre_spread,
                   climateRef = sim$climateRef_spread,
                   objFunCoresInternal = P(sim)$objFunCoresInternal,
                   tests = P(sim)$DEoptimTests,
                   maxFireSpread = P(sim)$maxFireSpread,
                   Nreps = P(sim)$objfunFireReps,
                   thresh = thresh,
                   .c = P(sim)$.c,
                   DEoptimControl = P(sim)$DEoptimControl,
                   .verbose = P(sim)$verbose,
                   visualizeDEoptim = P(sim)$visualizeDEoptim,
                   .plotSize = P(sim)$.plotSize,
                   .plots = P(sim)$.plots,
                   plotEvery = P(sim)$.plotInterval,
                   rep = P(sim)$.rep,
                   runName = runName,
                   progressFile = progressFile,
                   convergence = NULL, # set below: it depends on the start
                   sizeLik = P(sim)$sizeLik,
                   sizeLikDf = P(sim)$sizeLikDf,
                   weighted = P(sim)$weighted,
                   adWeight = P(sim)$adWeight,
                   link = spreadLink(P(sim)),
                   powerLawGamma = powerLawGammaOf(sim),
                   escapeSizeHa = escapeSizeHaOrNULL(P(sim)$escapeSizeHa),
                   jumpTries = P(sim)$jumpTries,
                   jumpMeanDist = P(sim)$jumpMeanDist,
                   yearAreaWeight = P(sim)$yearAreaWeight,
                   areaDistWeight = P(sim)$areaDistWeight,
                   penaliseRunaways = P(sim)$penaliseRunaways,
                   runawayEdgeFrac = P(sim)$runawayEdgeFrac,
                   runawayEdgeMin = P(sim)$runawayEdgeMin,
                   runawayBufferMultiple = runawayMultiple,
                   profileReps = if (diagnostics) P(sim)$profileReps else 0L,
                   simulateMembers = if (diagnostics) P(sim)$simulateMembers else 0L)
  lower <- P(sim)$lower
  n <- max(10L * length(lower), as.integer(P(sim)$nCoresNeeded))
  ## forced only when a start is drawn: a bad .ELFind is reported by the ledger's check, after the fit
  delayedAssign("seed", .elfSeed(sim$.ELFind))
  ## the population size runDEoptim() asks clusters for
  NP <- if (is.null(P(sim)$nCoresNeeded)) fireSenseUtils::deoptimPopulationSize(length(lower)) else
    as.integer(P(sim)$nCoresNeeded)
  emulated <- if (useEmulator(P(sim), length(lower)))
    Cache(emulatorStartPopulation, objArgs = emulatorObjectiveArgs(deArgs), lower = lower, upper = P(sim)$upper,
          NP = NP, control = P(sim)$emulatorControl, seed = seed,
          cores = best$cluster, nCoresNeeded = NP, libPath = deArgs$libPath, logPath = deArgs$logPath, runName = runName,
          .functionName = paste0("emulatorStartPopulation_", runName),
          ## the emulator's own code and the objective's, which Cache() does not digest
          .cacheExtra = list(emulatorBodies(), objectiveBodies()),
          omitArgs = c("cores", "libPath", "logPath", "runName"),
          useCache = P(sim)$useCache_DE)
  seedPop <- if (is.null(emulated)) P(sim)$initialpop else emulated$pop
  ## the early-stop rule in full (parameter, else options, else clusters' defaults; shorter from the emulator's
  ## start, convergenceOf()), so it is in the fit's key
  deArgs$convergence <- convergenceOf(P(sim)$convergence, seeded = !is.null(emulated))
  ## `initialpop` is NULL when there are no passing sets: DEoptim draws its start, and nothing in the cache
  ## key changes. With them it is a matrix (the passing sets, then uniform draws), a new fit with its own key,
  ## as is the retry below.
  fit <- function(initialpop = NULL) {
    args <- deArgs
    args$DEoptimControl <- deoptimControlStart(args$DEoptimControl, initialpop)
    Cache(do.call(runDEoptim, args),
        .functionName = fnName,
        ## The objective is a callee Cache() does not digest; its bodies are in the key, so a change in
        ## how fits are scored (e.g. fireSenseUtils' runaway censoring, 2026-10-01) is a cache miss,
        ## not an old fit served under the new rules.
        .cacheExtra = list(fnName, objectiveBodies()),
        ## runawayEdgeFrac/Min only change how quickly DEoptim moves away from an unlucky draw, not what a
        ## fit means, so fits cached under the 1-cell rule stay valid and must not rerun.
        omitArgs = c(".verbose", "cores", "paths", "logPath", "plotEvery", "progressFile", "runawayEdgeFrac", "runawayEdgeMin",
                     ## only the powerLaw link uses it
                     fireSenseUtils::powerLawOmitArgs(spreadLink(P(sim))),
                     ## NULL without an intercept (or with runawayBufferMultiple NA, or with the default spreadLink):
                     ## omitted then, so fits cached before it existed are found
                     fireSenseUtils::omitNullArgs(covCentre = sim$covCentre_spread,
                                                  climateRef = sim$climateRef_spread,
                                                  runawayBufferMultiple = runawayMultiple)),
        useCache = P(sim)$useCache_DE
    )
  }
  if (!is.null(seedPop)) passingPars <- NULL
  firstStart <- startPopulation(seedPop, passingPars, n, lower, P(sim)$upper, seed)
  ## clusters stops after the first generation when every member holds the fail value: DEoptim's
  ## uniform start had no member the objective accepts, and DE cannot improve on a population of fail
  ## values (ELF 6.6.2, 2026-10-20). Start again from a population the gates accept
  ## (feasibleInitialPopulation()); a second stop is final and says so.
  DE <- tryCatch(fit(initialpop = firstStart), deoptimNoFeasibleMember = function(e) {
    message("fireSense_spreadFit: ", conditionMessage(e),
            "\nDrawing a starting population that passes the spreadProb gates, as the threshold calibration does")
    gates <- spreadGates(covs, formulaToFit = sim$fireSense_spreadFormula, covMinMax = sim$covMinMax_spread,
                         covCentre = sim$covCentre_spread, climateRef = sim$climateRef_spread,
                         mutuallyExclusive = P(sim)$mutuallyExclusiveCols, maxFireSpread = P(sim)$maxFireSpread,
                         escapeSizeHa = escapeSizeHaOrNULL(P(sim)$escapeSizeHa), landscape = sim$rasterToMatch,
                         link = spreadLink(P(sim)), powerLawGamma = powerLawGammaOf(sim))
    localSingleThreadBLAS()
    pop <- feasibleInitialPopulation(n, lower = lower, upper = P(sim)$upper, gates = gates,
                                     seed = seed, cores = forkBudget(n))
    pop <- withPassingFirst(passingPars, pop, lower, P(sim)$upper)
    tryCatch(fit(initialpop = pop), deoptimNoFeasibleMember = function(e2)
      stop("fireSense_spreadFit: fit '", runName, "' failed again from a starting population of which ",
           attr(pop, "feasible"), " of ", nrow(pop), " members pass the spreadProb gates: ", conditionMessage(e2),
           call. = FALSE))
  })
  ## 1e6 is the objective's fail value (fireSenseUtils::.objfunSpreadFit). A population that is all
  ## fail values was never fitted; its "best" members are random draws.
  if (length(DE) && all(as.numeric(DE[[length(DE)]]$member$popval) >= 1e6))
    stop("fireSense_spreadFit: fit '", runName, "' failed: every member of the final DEoptim ",
         "population has the fail value 1e6, i.e. no parameter set passed the early-stop threshold (",
         thresh, " per fire) on the first two fire years.")
  DE
}

## The spread objective's bodies, for the Cache keys of the fit and the threshold calibration: the objective and
## everything it calls that decides its value (the spread link, which replicates are cut, the per-fire likelihood
## and score, and the terms added to it)
objectiveBodies <- function() {
  fns <- c("spreadProbFromIntegerCovs", "spreadProbCovsOfYear", "spreadProbFromCovs", "logisticAll", "powerLawLink",
           "dataYearBefore", "addBufferEdge", "bufferEdge", "bufferOpen", "sizeLikT", "kernelLogScore",
           "adStatistic", "adWeightAuto", "pooledFireSizes", "yearAreaNLL", "areaWeightedCvM",
           "pitMidRank", "adUniform")
  c(list(objfunSpreadFit = deparse(fireSenseUtils::.objfunSpreadFit),
         objFunInner = deparse(utils::getFromNamespace("objFunInner", "fireSenseUtils"))),
    stats::setNames(lapply(fns, function(f) deparse(utils::getFromNamespace(f, "fireSenseUtils"))), fns))
}

## The `link` the objective is given, from the parameters `Par` (`P(sim)`): "powerLaw" for `spreadLink = "powerLaw"`,
## else the `link` parameter's "logistic3pUpper", else NULL, its default, logistic3p
spreadLink <- function(Par) {
  if (identical(Par$spreadLink, "powerLaw")) "powerLaw" else if (identical(Par$link, "logistic3pUpper")) Par$link
}

## The bounds of `maxAsymptote`: for `spreadLink = "powerLaw"` it is p_c, which must stay above the floor and,
## per ELF, within `margin` of the calibrated threshold `p50` (fireSenseUtils::powerLawPcBounds())
maxAsymptoteBoundsOf <- function(spreadLink = "ceiling", p50 = NULL, margin = 0.005) {
  if (identical(spreadLink, "powerLaw") && !is.null(p50)) fireSenseUtils::powerLawPcBounds(p50, margin)
  else c(fireSenseUtils::spreadProbFloor, fireSenseUtils::spreadProbCeiling)
}

## Stop on a `spreadLink` the fit does not know, or one that cannot go with the other parameters
checkSpreadLink <- function(Par) {
  if (!Par$spreadLink %in% fireSenseUtils::spreadLinkOptions)
    stop("fireSense_spreadFit: spreadLink must be one of ", paste(dQuote(fireSenseUtils::spreadLinkOptions, FALSE), collapse = ", "),
         "; got ", dQuote(Par$spreadLink, FALSE))
  if (identical(Par$spreadLink, "powerLaw")) {
    if (identical(Par$link, "logistic3pUpper"))
      stop("fireSense_spreadFit: spreadLink = 'powerLaw' cannot be combined with link = 'logistic3pUpper'")
    if (!(is.numeric(Par$powerLawGamma) && length(Par$powerLawGamma) == 1L &&
          (is.na(Par$powerLawGamma) || (is.finite(Par$powerLawGamma) && Par$powerLawGamma > 0))))
      stop("fireSense_spreadFit: powerLawGamma must be one positive number, or NA to calibrate it; got ", format(Par$powerLawGamma))
  }
}

## `escapeSizeHa` as the objective takes it: NULL (the old fit) for NULL or NA
escapeSizeHaOrNULL <- function(x) if (length(x) && !is.na(x)) x

## `runawayBufferMultiple` as the objective takes it, by the same rule: NULL (a runaway scores the size it
## burned) for NULL or NA
runawayBufferMultipleOrNULL <- escapeSizeHaOrNULL

#' Turn the fit's re-score, profile and simulations into `sim` objects
#'
#' @param sim a `simList`.
#' @param DE the `runDEoptim()` result, with its attributes (not the reordered `sim$DE`).
#' @return the `simList`.
makeFitDiagnostics <- function(sim, DE) {
  if (is.null(DE)) return(sim)
  pop <- DE[[length(DE)]]$member$pop
  colnames(pop) <- names(P(sim)$lower)
  sim$spreadFitConvergence <- fireSenseUtils::fitConvergence(DE)

  scores <- attr(DE, "finalRescore")
  if (!is.null(scores)) {
    member <- seq_len(NROW(pop))
    reMean <- tapply(scores$value, scores$member, mean)
    reSD <- tapply(scores$value, scores$member, stats::sd)
    sim$spreadFitRescore <- data.table(member = member, pop,
                                       reMean = as.numeric(reMean[as.character(member)]),
                                       reSD = as.numeric(reSD[as.character(member)]))
  }

  ident <- fireSenseUtils::coefIdentifiability(pop, P(sim)$lower, P(sim)$upper)
  sim$spreadFitProfile <- attr(DE, "profile")
  if (!is.null(sim$spreadFitProfile)) {
    ident <- fireSenseUtils::identifiedInIsolation(ident, sim$spreadFitProfile)
    message("Identified in isolation (sign pinned by the population, and dropping it worsens the fit): ",
            paste(ident$coef[ident$identified], collapse = ", "),
            "\nNot identified in isolation: ", paste(ident$coef[!ident$identified], collapse = ", "))
  }
  sim$spreadFitIdentifiability <- ident

  sims <- attr(DE, "fitSims")
  if (!is.null(sims)) {
    sim$spreadFitSizes <- fireSenseUtils::scoreFireSizes(sims)
    sim$spreadFitLinkSaturation <- fireSenseUtils::linkSaturation(sims)
  }
  sim
}

## Year folds: every other year, in year order, so each fold spans the whole period
cvFolds <- function(years) {
  fold <- integer(length(years))
  fold[order(years)] <- rep_len(1:2, length(years))
  fold
}

#' Fit one cross-validation fold's complement and score its held-out years
#'
#' The per-fold work shared by [crossValidateSpread()] (both folds, one job) and
#' [crossValidateSpreadOneFold()] (one fold, one job -- the `heldOutFold` parameter). The fit's
#' `runName` (and so its `Cache` key, see `fitSpread()`) is suffixed by `k`, so fold 1, fold 2 and
#' a full fit (plain `P(sim)$.runName`, no suffix) never share a cache entry.
#'
#' The early-stop threshold is the full data's (`estimateFirstBlockThreshold()` on all years, the same
#' cached calibration as the full fit's). It is per fire, so it does not depend on how many fires the
#' two largest years of the fitted data hold. (It was a summed SNLL per year, which did: with the full
#' data's threshold, held-out fits of ELFs 4.3 and 5.2.1 (2026-09-29) never passed it and returned the
#' fail value for all 5000 generations, so each fold calibrated its own.)
#'
#' @param sim a `simList`.
#' @param covs `mod$covsX1000`.
#' @param fold integer vector, one element per year in `names(covs$historicalFires)`, from `cvFolds()`.
#' @param k integer; the fold (1 or 2) to hold out.
#' @return list: `sims`, this fold's simulated held-out years (`fireSenseUtils::simulateFireSizes()`), with
#'   `fold`; and `fit`, the fold's fit as a ledger row (`spreadFitLedgerRow()`) with all `simulateMembers` members.
fitAndScoreFold <- function(sim, covs, fold, k) {
  yearLists <- c("annualDTx1000", "fireBufferedListDT", "historicalFires")
  years <- names(covs$historicalFires)
  fitCovs <- heldCovs <- covs
  fitCovs[yearLists] <- lapply(covs[yearLists], function(x) x[years[fold != k]])
  heldCovs[yearLists] <- lapply(covs[yearLists], function(x) x[years[fold == k]])
  thresh <- estimateFirstBlockThreshold(sim, covs)$thresh
  DE <- fitSpread(sim, fitCovs, thresh, runName = paste0(P(sim)$.runName, "_cvFold", k), diagnostics = FALSE)
  best <- bestParamSets(DE, names(P(sim)$lower), n = max(1L, P(sim)$simulateMembers))
  s <- Cache(fireSenseUtils::simulateFireSizes, pop = best$params,
             fnArgs = spreadObjFunArgs(sim, heldCovs),
             .functionName = paste0("simulateHeldOut_", P(sim)$.runName, "_cvFold", k))
  ## the held-out years only: observed against simulated burning from the other fold's fit
  spreadFitValidationFigures(sim, as.matrix(best$params)[1L, ], heldCovs, paste0(P(sim)$.runName, "_cvFold", k, "_heldOut"))
  fit <- spreadFitLedgerRow(sim, length(DE), best$objFunVal, addFixedParsToLedger(best$params))
  list(sims = data.table(fold = k, s), fit = fit,
       fitYears = years[fold != k], heldOutYears = years[fold == k])
}

#' Two-fold cross-validation of the spread fit
#'
#' Fits the model to every other year and simulates the held-out years from the `simulateMembers`
#' best members of that fit, then the same the other way round. Uses the
#' same objective settings as the `run` event. Writes nothing to the ledger.
#'
#' @param sim a `simList`.
#' @param covs `mod$covsX1000`.
#' @return list: `sims` (from `fireSenseUtils::simulateFireSizes()`, with `fold`) and `score` (from
#'   `fireSenseUtils::scoreFireSizes()` on both folds together, so every year is predicted once).
crossValidateSpread <- function(sim, covs) {
  years <- names(covs$historicalFires)
  fold <- cvFolds(years)
  sims <- lapply(sort(unique(fold)), function(k) fitAndScoreFold(sim, covs, fold, k)$sims)
  sims <- rbindlist(sims)
  list(sims = sims, score = fireSenseUtils::scoreFireSizes(sims))
}

#' One cross-validation fold, run on its own (`heldOutFold`)
#'
#' Fits fold `k`'s complement and scores fold `k`'s held-out years only -- the other fold is never
#' run. Used by the `heldOutFold` parameter so the two folds can be separate jobs.
#'
#' @param sim a `simList`.
#' @param covs `mod$covsX1000`.
#' @param k integer; the fold (1 or 2) to hold out.
#' @return list: `sims` (this fold only), `score` (from `fireSenseUtils::scoreFireSizes()`), `fit` (see
#'   `fitAndScoreFold()`), `heldOutFold`, `fitYears` and `heldOutYears` (which fold, and its years), and the `formula` and `link` a prediction with `fit` needs.
crossValidateSpreadOneFold <- function(sim, covs, k) {
  years <- names(covs$historicalFires)
  fold <- cvFolds(years)
  res <- fitAndScoreFold(sim, covs, fold, k)
  list(sims = res$sims, score = fireSenseUtils::scoreFireSizes(res$sims), fit = res$fit,
       heldOutFold = k, fitYears = res$fitYears, heldOutYears = res$heldOutYears,
       formula = sim$fireSense_spreadFormula, link = P(sim)$link)
}

## The objective's arguments for simulating `covs`' years; the likelihood options do not matter, the
## escape rule does (the fit simulates every fire from its escape size)
spreadObjFunArgs <- function(sim, covs) {
  args <- list(landscape = sim$rasterToMatch,
       annualDTx1000 = covs$annualDTx1000, nonAnnualDTx1000 = covs$nonAnnualDTx1000,
       fireBufferedListDT = covs$fireBufferedListDT, historicalFires = covs$historicalFires,
       formulaToFit = sim$fireSense_spreadFormula, covMinMax = sim$covMinMax_spread,
       tests = P(sim)$DEoptimTests, maxFireSpread = P(sim)$maxFireSpread,
       objFunCoresInternal = P(sim)$objFunCoresInternal, Nreps = P(sim)$objfunFireReps,
       mutuallyExclusive = P(sim)$mutuallyExclusiveCols, doAssertions = FALSE, verbose = 0,
       link = spreadLink(P(sim)), escapeSizeHa = escapeSizeHaOrNULL(P(sim)$escapeSizeHa),
       jumpTries = P(sim)$jumpTries, jumpMeanDist = P(sim)$jumpMeanDist)
  ## the fit's centre and climate reference, not the held-out years': assigning NULL adds nothing, so without
  ## an intercept (or with the default spreadLink, or without powerLawGamma) the arguments (and the cache key of the simulation that
  ## holds them) are what they were
  args$covCentre <- sim$covCentre_spread
  args$climateRef <- sim$climateRef_spread
  args$powerLawGamma <- if (identical(args$link, "powerLaw")) powerLawGammaOf(sim)
  args
}
