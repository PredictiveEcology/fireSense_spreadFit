#' The threshold to fit with: Inf (no early stop) when calibration found none
#'
#' `pickThreshold()` returns `NA_real_` when every calibration trial failed. The threshold is only an
#' early-stop speed-up, so the correct fallback is no early stop, `Inf`. An `NA` left in place reaches
#' fireSenseUtils' objective, which compares `SNLL > thresh * 2` and dies on every cluster node with
#' "missing value where TRUE/FALSE needed" (ELF 6.2.1, heldOutFold 1). The conversion is done by the
#' consumer, not inside the cached calibration, so a cached `NA` is converted too.
#'
#' @param thresh the calibrated threshold: numeric, `NA`, or NULL (debug mode, passed through).
#' @param runName character; for the message.
#' @param nTrials integer; the number of calibration trials, for the message.
#' @return `thresh`, or `Inf` if it is `NA`.
noEarlyStopThreshold <- function(thresh, runName, nTrials) {
  if (!is.null(thresh) && length(thresh) == 1L && is.na(thresh)) {
    message("no SNLL threshold calibrated for ", runName, ": every one of ", nTrials,
            " trials failed; fitting without the early stop")
    return(Inf)
  }
  thresh
}
