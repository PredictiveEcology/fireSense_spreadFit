#' The calibrated early-stop threshold: `margin` times the best usable trial's first-block value per fire
#'
#' fireSenseUtils' objective bails out of a parameter set when the per-fire term of its first block of
#' years (the two largest fire years), divided by the number of fires in that block, is above `thresh`
#' (the objective's `firstBlockPerFire`). That is what `perFire` holds. Each fire's term is a kernel log
#' score that is 0 for a perfect prediction and never negative, so the threshold does not grow with the
#' number of fires, and `margin` times the best value is a margin above it.
#'
#' A trial is usable if it did not bail: it is finite, below `failVal`, and not `saturated` (its
#' spread probabilities failed the objective's "Too burny a landscape" / "Not spread out enough"
#' checks in the first block, so its fires scored the ceiling). The threshold is `margin` times the
#' smallest usable value. With `margin >= 1`, the best trial itself never bails under the result.
#'
#' @param perFire numeric; each trial's first-block value per fire, as run with no early stop.
#'   Non-numeric or non-finite entries (a crashed fork) count as failed.
#' @param saturated logical, recycled to `perFire`; the trial failed the first block's spread checks.
#' @param margin numeric; the `thresholdMargin` parameter.
#' @param failVal numeric; values at or above this are failed trials (fireSenseUtils' failVal).
#' @return the threshold, or `NA_real_` with a warning if no trial is usable
#'   (`noEarlyStopThreshold()` turns that into `Inf`).
pickThreshold <- function(perFire, saturated = FALSE, margin = 2, failVal = 1e6) {
  v <- suppressWarnings(as.numeric(perFire))
  ok <- is.finite(v) & v < failVal & !(rep_len(saturated, length(v)) %in% TRUE)
  if (!any(ok)) {
    warning("no threshold calibrated: no usable trial (", length(v), " tried)")
    return(NA_real_)
  }
  margin * min(v[ok])
}

#' One calibration trial: the objective with no early stop, reporting its first block
#'
#' Runs `fireSenseUtils::.objfunSpreadFit()` with `thresh = Inf` and `returnTerms = TRUE`, and reads
#' its `firstBlockPerFire` (the first block's per-fire term per fire, the quantity `thresh` is compared
#' with, see `pickThreshold()`) and `bailed` (a first-block year failed the spread checks).
#'
#' @param par,... passed to `.objfunSpreadFit()`; `...` must carry the fit's own objective settings,
#'   so the threshold is on the fit's scale.
#' @return `list(perFire =, saturated =)`; `perFire` is `NA` if the objective crashed, and then `error`
#'   holds the message.
trialFirstBlock <- function(par, ...) {
  res <- try(.objfunSpreadFit(par = par, thresh = Inf, returnTerms = TRUE, ...), silent = TRUE)
  if (inherits(res, "try-error"))
    return(list(perFire = NA_real_, saturated = NA, error = conditionMessage(attr(res, "condition"))))
  if (!"firstBlockPerFire" %in% names(res))
    return(list(perFire = NA_real_, saturated = NA, error = "the objective returned no firstBlockPerFire"))
  list(perFire = as.numeric(res[["firstBlockPerFire"]]), saturated = as.logical(res[["bailed"]]))
}

#' Stop if every calibration trial errored
#'
#' A threshold of `NA` (then `Inf`, see `noEarlyStopThreshold()`) is cached by `Cache()` and reused, so a
#' calibration whose every trial crashed (an objective that cannot run, not parameter sets that score
#' badly) must say so rather than return it.
#'
#' @param trials list of `trialFirstBlock()` results (an element that is not a list is a trial whose
#'   process died).
#' @return `NULL`, invisibly; stops with the first trial's error if all trials have one.
stopIfAllTrialsErrored <- function(trials) {
  errs <- vapply(trials, function(x) {
    if (!is.list(x)) "the trial's process died without a result"
    else if (is.null(x$error)) NA_character_ else x$error
  }, character(1))
  if (length(errs) && !anyNA(errs))
    stop("every calibration trial failed (", length(errs), " tried); the first error: ", errs[[1]], call. = FALSE)
  invisible(NULL)
}
