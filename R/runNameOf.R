#' The name of this fit
#'
#' Projects pass `.runName` as a parameter (`params$.globals`, as Biomass_core receives it);
#' older pipelines supply it as an input object. The parameter wins when it is a non-NA string,
#' otherwise the object `sim$.runName` is used, which may be `NULL`.
#'
#' @param sim A `simList`.
#' @return A character string, or `NULL` when neither the parameter nor the object is set.
runNameOf <- function(sim) {
  p <- P(sim)$.runName
  if (is.character(p) && length(p) == 1L && !is.na(p)) p else sim$.runName
}
