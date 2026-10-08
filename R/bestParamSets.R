#' The parameter sets the ledger records for this ELF
#'
#' The `n` distinct members of the final population with the lowest replicated mean
#' (`fireSenseUtils::runDEoptim()` attaches the re-score as `attr(DE, "finalRescore")`). Without a
#' re-score, the final population's own values are used. Not the `n` generations with the lowest best
#' value: those are copies of one, usually lucky, member.
#'
#' @param DE list of `DEoptim` objects returned by `fireSenseUtils::runDEoptim()`.
#' @param parNames character; parameter names, in the order of the columns of the population.
#' @param n integer; number of parameter sets to return.
#' @return list: `params` (matrix, one row per set) and `objFunVal` (their objective values), as
#'   returned by `fireSenseUtils::bestByReplicatedMean()`.
bestParamSets <- function(DE, parNames, n = 5L) {
  last <- DE[[length(DE)]]$member
  pop <- last$pop
  colnames(pop) <- parNames
  scores <- attr(DE, "finalRescore")
  if (is.null(scores))
    scores <- data.table::data.table(member = seq_len(NROW(pop)), rep = 1L, value = as.numeric(last$popval))
  best <- fireSenseUtils::bestByReplicatedMean(pop, scores, n = n)
  list(params = best$params, objFunVal = best$objFunVal)
}
