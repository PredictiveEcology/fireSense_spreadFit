#' The ledger row of a fit
#'
#' One row of the shared fit ledger: the polygon's geometry, the columns in
#' `fireSenseUtils::spreadFitAdditionalColNamesTxt`, and `polygonID`. The `run` event writes it to the
#' ledger, and a held-out fold saves it beside its held-out result, so `fireSense_spreadPredict` reads
#' either the same way.
#'
#' @param sim a `simList`; uses `.ELFind`, `studyArea`, `sppEquiv`, `nonForestedLCCGroups`,
#'   `missingLCCgroup` and `covMinMax_spread`.
#' @param numIterations integer, the number of DEoptim generations.
#' @param objFunVal numeric, the objective values of the members in `params`.
#' @param params `data.table`, one row per member, with `hillSlope1` and `inflectionPoint1` (see `addFixedParsToLedger()`).
#' @return an `sf` object with one row.
spreadFitLedgerRow <- function(sim, numIterations, objFunVal, params) {
  ## covMinMax_spread: prediction rescales covariates with it, exactly as this fit did
  df <- data.frame(I(list(numIterations)),
                   I(list(objFunVal)),
                   I(list(params)),
                   I(list(sim$sppEquiv)),
                   I(list(sim$nonForestedLCCGroups)),
                   I(list(sim$missingLCCgroup)),
                   I(list(sim$covMinMax_spread))) |>
    setNames(fireSenseUtils::spreadFitAdditionalColNamesTxt)
  # The ledger is keyed by polygon identity, NOT by run label -- see the
  # `.ELFind` input declaration. This row is shared cloud state that every
  # other project reads, so validate before writing.
  polygonID <- sim$.ELFind
  if (!is.character(polygonID) || length(polygonID) != 1L ||
      is.na(polygonID) || !nzchar(polygonID))
    stop("fireSense_spreadFit: `sim$.ELFind` must be a single non-empty character ",
         "identifying the polygon being fit; got: ",
         paste(format(polygonID), collapse = ", "))
  df <- data.frame(df, "polygonID" = polygonID)

  saHere <- if (is(sim$studyArea, "SpatVector")) sf::st_as_sf(sim$studyArea) else sim$studyArea
  saHere <- sf::st_as_sf(sf::st_geometry(saHere))
  sf::st_geometry(saHere) <- "geometry"
  dplyr::mutate(saHere, df)
}
