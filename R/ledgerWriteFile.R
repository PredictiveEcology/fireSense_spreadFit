#' The ledger file a fit is written to
#'
#' `"latest"` is how readers find fits, not a file: a fit goes to the file named for its own fire
#' years and model (`fireSenseUtils::spreadFitFilenameFor()`), which is then the latest one that has
#' this polygon. Any other `spreadFitFilename` is used as it is.
#'
#' @param spreadFitFilename The module's `spreadFitFilename` parameter.
#' @param fireYears fireSense_dataPrepFit's `fireYears`, the window the fit used; `NULL` when that
#'   module is not in the simulation.
#' @param annualNames Names of `sim$fireSense_annualSpreadFitCovariates` (`"year2001"`, ...), for
#'   the window when `fireYears` is not available.
#' @return A file name.
ledgerWriteFile <- function(spreadFitFilename, fireYears, annualNames) {
  if (!identical(spreadFitFilename, "latest"))
    return(spreadFitFilename)
  if (!length(fireYears) || all(is.na(fireYears)))
    fireYears <- suppressWarnings(as.integer(gsub("\\D", "", annualNames)))
  fireSenseUtils::spreadFitFilenameFor(fireYears)
}
