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

#' The ledger file this fit writes
#'
#' `ledgerWriteFile()` for a run: the module's `spreadFitFilename`, the fire years from
#' `fireSense_dataPrepFit`, and the window of the annual covariates.
#'
#' @param sim A `simList`.
#' @return A file name.
spreadFitLedgerFile <- function(sim) {
  ledgerWriteFile(P(sim)$spreadFitFilename,
                  P(sim, module = "fireSense_dataPrepFit")$fireYears,
                  names(sim$fireSense_annualSpreadFitCovariates))
}

#' The folder this fit's figures, held-out results and DEoptim plots go to
#'
#' The parameter `fitOutputPath` when set; otherwise next to this fit's ledger file, a folder named
#' for the polygon (`sim$.ELFind`) and that file (`fireSenseUtils::fitOutputPath()`). A fit depends
#' only on those, so its outputs do not go with whichever scenario and replicate ran it.
#'
#' @param sim A `simList`.
#' @return A path. Nothing is created.
spreadFitOutputPath <- function(sim) {
  if (!is.null(P(sim)$fitOutputPath)) return(P(sim)$fitOutputPath)
  fireSenseUtils::fitOutputPath(inputPath(sim), sim$.ELFind, spreadFitLedgerFile(sim))
}
