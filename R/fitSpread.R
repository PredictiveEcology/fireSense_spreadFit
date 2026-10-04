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
#' @param runName character; labels the run and its cache entry.
#' @param diagnostics logical; `FALSE` skips the profile and the simulations.
#' @return the `runDEoptim()` result. Stops if every member of the final population has the
#'   objective's fail value: no parameter set ever passed `thresh`, so there is no fit.
fitSpread <- function(sim, covs, thresh, runName, diagnostics = TRUE) {
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
  DE <- Cache(runDEoptim(landscape = sim$rasterToMatch,
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
                   sizeLik = P(sim)$sizeLik,
                   sizeLikDf = P(sim)$sizeLikDf,
                   weighted = P(sim)$weighted,
                   adWeight = P(sim)$adWeight,
                   link = spreadLink(P(sim)$link),
                   escapeSizeHa = escapeSizeHaOrNULL(P(sim)$escapeSizeHa),
                   jumpTries = P(sim)$jumpTries,
                   jumpMeanDist = P(sim)$jumpMeanDist,
                   yearAreaWeight = P(sim)$yearAreaWeight,
                   areaDistWeight = P(sim)$areaDistWeight,
                   penaliseRunaways = P(sim)$penaliseRunaways,
                   runawayEdgeFrac = P(sim)$runawayEdgeFrac,
                   runawayEdgeMin = P(sim)$runawayEdgeMin,
                   profileReps = if (diagnostics) P(sim)$profileReps else 0L,
                   simulateMembers = if (diagnostics) P(sim)$simulateMembers else 0L),
        .functionName = fnName,
        .cacheExtra = fnName,
        ## runawayEdgeFrac/Min only change how quickly DEoptim moves away from an unlucky draw, not what a
        ## fit means, so fits cached under the 1-cell rule stay valid and must not rerun.
        omitArgs = c(".verbose", "cores", "paths", "logPath", "plotEvery", "runawayEdgeFrac", "runawayEdgeMin"),
        useCache = P(sim)$useCache_DE
  )
  ## 1e6 is the objective's fail value (fireSenseUtils::.objfunSpreadFit). A population that is all
  ## fail values was never fitted; its "best" members are random draws.
  if (length(DE) && all(as.numeric(DE[[length(DE)]]$member$popval) >= 1e6))
    stop("fireSense_spreadFit: fit '", runName, "' failed: every member of the final DEoptim ",
         "population has the fail value 1e6, i.e. no parameter set passed the SNLL threshold (",
         thresh, ") on the first two fire years.")
  DE
}

## The `link` the objective is given: NULL is its default, logistic3p
spreadLink <- function(link) if (identical(link, "logistic3pUpper")) link

## `escapeSizeHa` as the objective takes it: NULL (the old fit) for NULL or NA
escapeSizeHaOrNULL <- function(x) if (length(x) && !is.na(x)) x

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
#' a full fit (plain `sim$.runName`, no suffix) never share a cache entry.
#'
#' The SNLL threshold is calibrated on the fitted years (`estimateSNLLThresholdPostLargeFires()`), not
#' taken from the full data: it bounds the SNLL of the two largest fire years, and a fold's two largest
#' years are not the full data's. With the full data's threshold, held-out fits of ELFs 4.3 and 5.2.1
#' (2026-09-29) never passed it and returned the fail value for all 5000 generations.
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
  thresh <- estimateSNLLThresholdPostLargeFires(sim, fitCovs)
  DE <- fitSpread(sim, fitCovs, thresh, runName = paste0(sim$.runName, "_cvFold", k), diagnostics = FALSE)
  best <- bestParamSets(DE, names(P(sim)$lower), n = max(1L, P(sim)$simulateMembers))
  s <- Cache(fireSenseUtils::simulateFireSizes, pop = best$params,
             fnArgs = spreadObjFunArgs(sim, heldCovs),
             .functionName = paste0("simulateHeldOut_", sim$.runName, "_cvFold", k))
  ## the held-out years only: observed against simulated burning from the other fold's fit
  spreadFitValidationFigures(sim, as.matrix(best$params)[1L, ], heldCovs, paste0(sim$.runName, "_cvFold", k, "_heldOut"))
  fit <- spreadFitLedgerRow(sim, length(DE), best$objFunVal, addHillSlope1ToLedger(best$params))
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
  list(landscape = sim$rasterToMatch,
       annualDTx1000 = covs$annualDTx1000, nonAnnualDTx1000 = covs$nonAnnualDTx1000,
       fireBufferedListDT = covs$fireBufferedListDT, historicalFires = covs$historicalFires,
       formulaToFit = sim$fireSense_spreadFormula, covMinMax = sim$covMinMax_spread,
       tests = P(sim)$DEoptimTests, maxFireSpread = P(sim)$maxFireSpread,
       objFunCoresInternal = P(sim)$objFunCoresInternal, Nreps = P(sim)$objfunFireReps,
       mutuallyExclusive = P(sim)$mutuallyExclusiveCols, doAssertions = FALSE, verbose = 0,
       link = spreadLink(P(sim)$link), escapeSizeHa = escapeSizeHaOrNULL(P(sim)$escapeSizeHa),
       jumpTries = P(sim)$jumpTries, jumpMeanDist = P(sim)$jumpMeanDist)
}
