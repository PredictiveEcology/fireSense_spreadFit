#' The final DEoptim population, as a table
#'
#' One row per member of the last generation: the fitted parameters (named as `parNames`, in the order of the
#' population's columns) and `objFunVal`, the objective value DEoptim holds for it.
#'
#' @param DE list of `DEoptim` objects returned by `fireSenseUtils::runDEoptim()`.
#' @param parNames character; parameter names, as for [bestParamSets()].
#' @return a `data.table`.
finalPopulation <- function(DE, parNames) {
  last <- DE[[length(DE)]]$member
  pop <- data.table::as.data.table(last$pop)
  data.table::setnames(pop, parNames)
  pop[, objFunVal := as.numeric(last$popval)]
  pop[]
}

#' Write the final population of a fit
#'
#' To `DEoptimPopulation_<runName>.csv` in the module's output folder, beside the fit's progress file. The
#' population is the pool [scoreTermsPool()] (mode "terms") reads.
#'
#' @inheritParams finalPopulation
#' @param sim a `simList`.
#' @param runName character; names the file.
#' @return the file's path, invisibly.
writeFinalPopulation <- function(sim, DE, parNames, runName) {
  out <- file.path(checkPath(file.path(outputPath(sim), currentModule(sim)), create = TRUE),
                   paste0("DEoptimPopulation_", runName, ".csv"))
  data.table::fwrite(finalPopulation(DE, parNames), out)
  message("fireSense_spreadFit: wrote the final population to ", out)
  invisible(out)
}

#' The parameter vectors mode "terms" scores
#'
#' `termsPool` is `NA` (default): the sets of this ELF's ledger row ([ledgerFitRow()]). Otherwise a csv with a column
#' for every fitted parameter, e.g. `DEoptimPopulation_<runName>.csv` ([writeFinalPopulation()]) or any table of
#' vectors; other columns are ignored.
#'
#' @inheritParams ledgerFitRow
#' @return a `data.table` of the fitted parameters, one row per vector.
termsPoolSets <- function(sim) {
  pool <- P(sim)$termsPool
  if (isTRUE(is.na(pool))) {
    fit <- ledgerFitRow(sim, "terms")
    return(fit$sets[, fit$parNames, with = FALSE])
  }
  parNames <- names(P(sim)$lower)
  sets <- data.table::fread(pool)
  if (!all(parNames %in% names(sets)))
    stop("fireSense_spreadFit: `termsPool` ", pool, " does not hold this run's parameters: ",
         paste(setdiff(parNames, names(sets)), collapse = ", "), " missing", call. = FALSE)
  sets[, parNames, with = FALSE]
}

#' Mode "terms": the objective's separate terms for a pool of parameter vectors
#'
#' Every vector of [termsPoolSets()] is evaluated `termsReps` times on the fit's data and objective
#' (`spreadFitDEArgs()`, on a cluster built as the fit's, `withObjectiveCluster()`), each with its terms returned
#' separately (`fireSenseUtils::scorePopulationTerms()`), and replicate `r` of every vector uses the same random
#' numbers (`termsSeed + r`), so vectors differ by their parameters, not by luck. The terms can then be weighted
#' afterwards without simulating again. Nothing is fitted and the ledger is not written.
#'
#' @param sim a `simList`, after `spreadFitPrepare`.
#' @param covs `mod$covsX1000`.
#' @param runName character; labels the cluster and the output file.
#' @return a `data.table`: `member` (row of the pool), `rep`, the parameters and the terms. Also written to
#'   `file.path(outputPath(sim), currentModule(sim), "spreadFitTerms_<runName>.csv")`.
scoreTermsPool <- function(sim, covs, runName) {
  sets <- termsPoolSets(sim)
  deArgs <- spreadFitDEArgs(sim, covs, runName, diagnostics = FALSE)
  objArgs <- emulatorObjectiveArgs(deArgs)
  message("fireSense_spreadFit: scoring ", nrow(sets), " parameter vectors, ", P(sim)$termsReps, " reps each")
  terms <- withObjectiveCluster(objArgs, cores = deArgs$cores, nCoresNeeded = spreadFitNP(sim), libPath = deArgs$libPath,
                                logPath = deArgs$logPath, runName = paste0("terms_", runName), FUN = function(cl)
    fireSenseUtils::scorePopulationTerms(sets, fn = fireSenseUtils::.objfunSpreadFit, reps = as.integer(P(sim)$termsReps),
                                         cl = cl, seed = as.integer(P(sim)$termsSeed),
                                         fnArgs = objectiveFnArgs(objArgs, cl)))
  out <- cbind(terms[, c("member", "rep")], sets[terms$member], terms[, -c("member", "rep")])
  file <- file.path(checkPath(file.path(outputPath(sim), currentModule(sim)), create = TRUE),
                    paste0("spreadFitTerms_", runName, ".csv"))
  data.table::fwrite(out, file)
  message("fireSense_spreadFit: wrote the terms to ", file)
  out
}
