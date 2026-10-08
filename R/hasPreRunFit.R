#' Does the shared ledger already hold a fit for the polygon this run fits?
#'
#' `sim$studyAreaWithSpreadParams` is whatever ledger rows intersect the study
#' area: neighbours' fits, this polygon's fit, or nothing at all. Its class
#' therefore says nothing about *this* polygon. Deciding by class (`is(x, "sf")`)
#' skipped the fit for every unfitted ELF -- `fireSense_dataPrepFit` assigns the
#' object whenever the ledger loads -- and left `fireSense_spreadPredict` to die on
#' an empty parameter set (2026-09-07, every job of the 41-ELF batch).
#'
#' @param sim A `simList`. Uses `sim$.ELFind` as the polygon identity (the ledger
#'   key), falling back to `sim$.runName`.
#' @return Logical, length one.
hasPreRunFitForThisPolygon <- function(sim) {
  sa <- sim$studyAreaWithSpreadParams
  if (!(is(sa, "sf") || is(sa, "data.frame")) || NROW(sa) == 0L) return(FALSE)
  if (!"polygonID" %in% names(sa)) return(FALSE)
  id <- sim$.ELFind
  if (is.null(id) || !length(id)) id <- sim$.runName
  if (is.null(id) || !length(id)) return(FALSE)
  isTRUE(any(as.character(sa$polygonID) == as.character(id[1])))
}
