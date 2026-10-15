#' Does this run fit the polygon, or use the ledger's fit? Stops when `doFit` is FALSE and there is no fit.
#'
#' Reads the module's `doFit` and its deprecated aliases through `fireSenseUtils::resolveDoFit()` and
#' decides with `fireSenseUtils::mustFit()`.
#'
#' @param sim a `simList`.
#' @param message logical; send the deprecation message (from `init`, once).
#' @return Logical, length one: TRUE to fit, FALSE to use the ledger's fit.
mustFitThisPolygon <- function(sim, message = FALSE) {
  doFit <- fireSenseUtils::resolveDoFit(P(sim)$doFit, P(sim)$refitExisting, P(sim)$stopIfNoPreRunFit,
                                        module = "fireSense_spreadFit", message = message)
  fireSenseUtils::mustFit(doFit, hasPreRunFitForThisPolygon(sim), sim$.ELFind)
}
