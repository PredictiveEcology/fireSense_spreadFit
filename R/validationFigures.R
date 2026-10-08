#' Figures comparing a spread fit with its data
#'
#' Simulates the fire years in `covs` from `par`, once, with the objective's settings and
#' `objfunFireReps` replicates (`fireSenseUtils::spreadFitValidationData()`), and draws two
#' figures from that one table with `Plots()`:
#'
#' * `spreadFitObservedVsSimulated_<label>`: the observed and the simulated burned share of
#'   pixel-years, binned by each covariate (`fireSenseUtils::plotSpreadFitValidation()`). This is
#'   the figure that can show misfit.
#' * `spreadFitResponseCurves_<label>`: the fitted response curves
#'   (`fireSenseUtils::plotSpreadFitResponse()`), the model's response rather than a fit to data.
#'
#' Nothing is simulated unless `.plots` asks for a plot. The simulation costs about one objective
#' evaluation (no early stop): 18 s on ELF 5.3.2 (36 fire years, 1.0 million pixel-years, 50
#' replicates).
#'
#' @param sim a `simList`.
#' @param par named numeric; the parameter set, named as `P(sim)$lower`.
#' @param covs the covariates as integers x 1000, restricted to the years to show.
#' @param label character; ends the file names.
#' @return `NULL`, invisibly.
spreadFitValidationFigures <- function(sim, par, covs, label) {
  if (!anyPlotting(P(sim)$.plots)) return(invisible(NULL))
  d <- do.call(spreadFitValidationData,
               c(list(par = par, seed = .elfSeed(sim$.ELFind)), spreadObjFunArgs(sim, covs)))
  Plots(d, fn = plotSpreadFitValidation, types = P(sim)$.plots, path = figurePath(sim),
        filename = paste0("spreadFitObservedVsSimulated_", label),
        ggsaveArgs = list(width = 12, height = 9))
  Plots(d, fn = plotSpreadFitResponse, types = P(sim)$.plots, path = figurePath(sim),
        filename = paste0("spreadFitResponseCurves_", label),
        ggsaveArgs = list(width = 12, height = 9))
  invisible(NULL)
}
