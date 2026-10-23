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

#' Does this job use the `"powerLaw"` calibration?
#'
#' Only a job that fits (`mustFitThisPolygon()`), validates (`mode` has "validate") or runs a held-out fold
#' (`heldOutFold`) sets bounds or simulates fires with the calibrated `p_c50` and `gamma`. A job that uses the
#' ledger's fit reads them from the ledger row, so it skips the calibration.
#'
#' @param sim a `simList`.
#' @return Logical, length one.
needsPowerLawCalibration <- function(sim) {
  heldOut <- !isTRUE(is.na(P(sim)$heldOutFold))
  heldOut || "validate" %in% P(sim)$mode || mustFitThisPolygon(sim)
}
