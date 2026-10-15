#' Should this run fit the polygon, use the ledger's fit, or stop?
#'
#' `doFit` is `NA` (fit only when the ledger has no fit for this polygon), `TRUE` (always fit,
#' replacing a ledger fit) or `FALSE` (never fit: use the ledger fit, or stop if there is none).
#' The deprecated `refitExisting` (TRUE is `doFit = TRUE`) and `stopIfNoPreRunFit` (TRUE is
#' `doFit = FALSE`; FALSE only allowed a fit, which is what `NA` does) are mapped here, so
#' `init` and `run` read one value. An explicit `doFit` is not overridden by an alias.
#'
#' @param par the module's parameters, `P(sim)`.
#' @param message logical; send the deprecation message (once, from `init`).
#' @return `NA`, `TRUE` or `FALSE`.
resolveDoFit <- function(par, message = FALSE) {
  aliasSet <- !is.na(par$refitExisting) || !is.na(par$stopIfNoPreRunFit)
  if (aliasSet && message)
    message("fireSense_spreadFit: parameters `refitExisting` and `stopIfNoPreRunFit` are deprecated; ",
            "use `doFit` (NA: fit only if the ledger has no fit; TRUE: always fit; FALSE: never fit).")
  if (!is.na(par$doFit)) return(isTRUE(par$doFit))
  if (isTRUE(par$refitExisting)) return(TRUE)
  if (isTRUE(par$stopIfNoPreRunFit)) return(FALSE)
  NA
}

#' Does this run fit the polygon? Stops when `doFit` forbids the fit the ledger needs.
#'
#' @param sim a `simList`.
#' @param message logical; passed to [resolveDoFit()].
#' @return Logical, length one: TRUE to fit, FALSE to use the ledger's fit.
mustFitThisPolygon <- function(sim, message = FALSE) {
  doFit <- resolveDoFit(P(sim), message = message)
  if (isTRUE(doFit)) return(TRUE)
  hasFit <- hasPreRunFitForThisPolygon(sim)
  if (!hasFit && isFALSE(doFit))
    stop("There is no pre-run SpreadFit (sim$studyAreaWithSpreadParams) for polygon ", sim$.ELFind[1],
         " and `doFit` is FALSE. Set `doFit = NA` (or TRUE) to fit it.")
  !hasFit
}
