#' The calibrated SNLL threshold: `margin` times the best usable trial's first-block SNLL
#'
#' fireSenseUtils' objective bails out of a parameter set when the SNLL of its first block of years
#' (the two largest fire years) is above `thresh * <years in the block>`, i.e. when the block's
#' average annual SNLL is above `thresh`. The unit of `thresh` is therefore the objective's own
#' average annual SNLL of the first block (the objective's `firstBlockSNLL`), and that is what `annual` holds.
#'
#' A trial is usable if it did not bail: it is finite, below `failVal`, and not `saturated` (its
#' spread probabilities failed the objective's "Too burny a landscape" / "Not spread out enough"
#' checks in the first block, so it scored the minLik floor). The threshold is `margin` times the
#' smallest usable value, rounded up because the objective rounds `thresh`. With `margin >= 1` and a
#' positive best value, the best trial itself never bails under the result.
#'
#' @param annual numeric; each trial's first-block average annual SNLL, as run with no early stop.
#'   Non-numeric or non-finite entries (a crashed fork) count as failed.
#' @param saturated logical, recycled to `annual`; the trial failed the first block's spread checks.
#' @param margin numeric; the `thresholdMargin` parameter.
#' @param failVal numeric; values at or above this are failed trials (fireSenseUtils' failVal).
#' @return the threshold, or `NA_real_` with a warning if no trial is usable
#'   (`noEarlyStopThreshold()` turns that into `Inf`).
pickThreshold <- function(annual, saturated = FALSE, margin = 2, failVal = 1e6) {
  v <- suppressWarnings(as.numeric(annual))
  ok <- is.finite(v) & v < failVal & !(rep_len(saturated, length(v)) %in% TRUE)
  if (!any(ok)) {
    warning("no threshold calibrated: no usable trial (", length(v), " tried)")
    return(NA_real_)
  }
  ceiling(margin * min(v[ok]))
}

#' One calibration trial: the objective with no early stop, reporting its first block
#'
#' Runs `fireSenseUtils::.objfunSpreadFit()` with `thresh = Inf` and `returnTerms = TRUE`, and reads
#' its `firstBlockSNLL` (the first block's average annual SNLL, the quantity `thresh` is compared with,
#' see `pickThreshold()`) and `bailed` (a first-block year failed the spread checks).
#'
#' @param par,... passed to `.objfunSpreadFit()`; `...` must carry the fit's own objective settings,
#'   so the threshold is on the fit's scale.
#' @return `list(annual =, saturated =)`; `annual` is `NA` if the objective crashed.
trialFirstBlock <- function(par, ...) {
  res <- try(.objfunSpreadFit(par = par, thresh = Inf, returnTerms = TRUE, ...), silent = TRUE)
  if (inherits(res, "try-error") || !"firstBlockSNLL" %in% names(res))
    return(list(annual = NA_real_, saturated = NA))
  list(annual = as.numeric(res[["firstBlockSNLL"]]), saturated = as.logical(res[["bailed"]]))
}
