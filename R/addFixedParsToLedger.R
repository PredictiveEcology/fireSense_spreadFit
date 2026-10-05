#' Insert the fixed `hillSlope1` and `inflectionPoint1` into the ledger's parameter sets
#'
#' `hillSlope1` (the spread link's slope) and `inflectionPoint1` (its Richards exponent) are fixed at
#' 1, not fitted (`fireSenseUtils::fixedLogisticPars`; see `estimateSpreadParams()`): neither is
#' identifiable together with the covariate coefficients, so they are no longer part of `DEoptim`'s
#' parameter space, and `paramsBest` (from `bestParamSets()`) does not include them either.
#' `fireSense_spreadPredict` splits a ledger row's parameters from its covariates BY NAME, and the
#' spread link (`fireSenseUtils::logistic3p()`/`logistic3pUpper()`) reads the result BY POSITION --
#' `maxAsymptote`, `hillSlope1`, `inflectionPoint1`, and (with the upper-tail link) `upperTail1`.
#' This restores those positions for new fits, so a new ledger row predicts with both at 1 exactly
#' like an old row whose `hillSlope1` and `inflectionPoint1` happened to be fitted at that value, and
#' an old row keeps predicting with its own fitted values.
#'
#' @param paramsBest a `data.table`, one row per parameter set, as `bestParamSets()$params` returns
#'   it: columns named `names(P(sim)$lower)`, so without `hillSlope1` and `inflectionPoint1`,
#'   `maxAsymptote` first.
#' @return `paramsBest` with `hillSlope1` and `inflectionPoint1` columns of 1s inserted right after
#'   `maxAsymptote`.
addFixedParsToLedger <- function(paramsBest) {
  paramsBest <- data.table::copy(paramsBest)
  fixedPars <- fireSenseUtils::fixedLogisticPars
  for (nm in names(fixedPars)) data.table::set(paramsBest, NULL, nm, fixedPars[[nm]])
  data.table::setcolorder(paramsBest,
                          append(setdiff(names(paramsBest), names(fixedPars)), names(fixedPars), after = 1L))
  paramsBest
}
