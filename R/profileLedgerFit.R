#' Profile the ledger's fit of this ELF, without fitting
#'
#' Mode `"profile"`: the objective is evaluated at the best parameter set of this ELF's ledger row (its first
#' set, the best by replicated mean), and at that set with each coefficient of `profileCoefs` in turn set to each
#' value of `profileAt`, `profileReps` times each on paired seeds (`fireSenseUtils::profileCoefficients()`). The
#' values may lie outside the fit's bounds, e.g. below the fuel coefficients' lower bound of 0. The objective and
#' its data are the fit's (`spreadFitDEArgs()`, `emulatorObjectiveArgs()`), on a cluster built as the fit's
#' (`withObjectiveCluster()`). Nothing is fitted and the ledger is not written.
#'
#' Stops when the ledger has no row for this ELF, or when the row was fitted on other parameters or other data
#' than this run would fit: other parameter names, another covariate centre, or another power-law `gamma`. A
#' profile is only about the fit when the objective is the fit's.
#'
#' @param sim a `simList`, after `spreadFitPrepare`.
#' @param covs `mod$covsX1000`.
#' @param runName character; labels the cluster and the output file.
#' @return the `profileCoefficients()` table, with `ledgerObjFunVal` (the ledger's replicated mean of the best set,
#'   to compare with the reference row's `mean`). It is also written to
#'   `file.path(outputPath(sim), currentModule(sim), "spreadFitProfile_<runName><fitFileSuffix>.csv")`.
profileLedgerFit <- function(sim, covs, runName) {
  if (!hasPreRunFitForThisPolygon(sim))
    stop("fireSense_spreadFit: mode 'profile' profiles this ELF's ledger fit, but the ledger has no row for ",
         paste(sim$.ELFind, collapse = ", "), call. = FALSE)
  if (all(is.na(P(sim)$profileAt)))
    stop("fireSense_spreadFit: mode 'profile' needs `profileAt`, the values to profile each coefficient at",
         call. = FALSE)
  sa <- sim$studyAreaWithSpreadParams
  row <- sa[as.character(sa$polygonID) %in% as.character(sim$.ELFind), ][1, ]
  parNames <- names(P(sim)$lower)
  ledgerPars <- as.data.frame(row$params[[1]])
  if (!all(parNames %in% names(ledgerPars)))
    stop("fireSense_spreadFit: the ledger row's parameters (", paste(names(ledgerPars), collapse = ", "),
         ") are not this run's (", paste(parNames, collapse = ", "), "); it was fitted with another formula",
         call. = FALSE)
  best <- unlist(ledgerPars[1L, parNames])
  ## the same data and link as the fit, or the profile is not about the fit
  same <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-6))
  centre <- sim$covCentre_spread
  rowCentre <- row$covCentre_spread[[1]]
  if (!is.null(centre) && !same(unlist(rowCentre[names(centre)]), unlist(centre)))
    stop("fireSense_spreadFit: this run's covariate centre differs from the ledger fit's, so its data differ",
         call. = FALSE)
  if (P(sim)$spreadLink %in% fireSenseUtils::spreadPowerLawLinks &&
      fireSenseUtils::spreadFitPowerLawGammaTxt %in% names(row) &&
      !same(row[[fireSenseUtils::spreadFitPowerLawGammaTxt]], powerLawGammaOf(sim)))
    stop("fireSense_spreadFit: this run's power-law gamma (", powerLawGammaOf(sim), ") differs from the ledger fit's (",
         row[[fireSenseUtils::spreadFitPowerLawGammaTxt]], ")", call. = FALSE)

  covCoefs <- fireSenseUtils::spreadDesignCols(sim$fireSense_spreadFormula) # the intercept, if any, and the terms
  coefs <- if (all(is.na(P(sim)$profileCoefs))) covCoefs else
    unique(unlist(lapply(P(sim)$profileCoefs, function(x)
      if (identical(x, "fuel")) intersect(covCoefs, fuelColumns(sim$fireSense_nonAnnualSpreadFitCovariates)) else x)))
  if (!length(coefs) || !all(coefs %in% covCoefs))
    stop("fireSense_spreadFit: `profileCoefs` must name covariate coefficients of the fit (", paste(covCoefs, collapse = ", "),
         ") or be \"fuel\"; got ", paste(P(sim)$profileCoefs, collapse = ", "), call. = FALSE)
  at <- stats::setNames(rep(list(sort(unique(as.numeric(P(sim)$profileAt)))), length(coefs)), coefs)

  deArgs <- spreadFitDEArgs(sim, covs, runName, diagnostics = FALSE)
  objArgs <- emulatorObjectiveArgs(deArgs)
  NP <- if (is.null(P(sim)$nCoresNeeded)) fireSenseUtils::deoptimPopulationSize(length(parNames)) else
    as.integer(P(sim)$nCoresNeeded)
  message("fireSense_spreadFit: profiling the ledger fit of ", paste(sim$.ELFind, collapse = ", "), ": ",
          paste(coefs, collapse = ", "), " at ", paste(at[[1]], collapse = ", "), ", ", P(sim)$profileReps, " reps each")
  prof <- withObjectiveCluster(objArgs, cores = deArgs$cores, nCoresNeeded = NP, libPath = deArgs$libPath,
                               logPath = deArgs$logPath, runName = paste0("profile_", runName), FUN = function(cl)
    fireSenseUtils::profileCoefficients(best, pop = NULL, fn = fireSenseUtils::.objfunSpreadFit,
                                        reps = as.integer(P(sim)$profileReps), cl = cl,
                                        fnArgs = if (is.null(cl)) objArgs else objArgs[setdiff(names(objArgs), emulatorDataArgs)],
                                        coefs = coefs, at = at))
  prof[, ledgerObjFunVal := as.numeric(row$objFunVal[[1]][1])]
  out <- file.path(checkPath(file.path(outputPath(sim), currentModule(sim)), create = TRUE),
                   paste0("spreadFitProfile_", runName, P(sim)$fitFileSuffix, ".csv"))
  data.table::fwrite(prof, out)
  message("fireSense_spreadFit: wrote the profile to ", out)
  prof
}
