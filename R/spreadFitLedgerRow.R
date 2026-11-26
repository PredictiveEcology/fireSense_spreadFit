#' The ledger row of a fit
#'
#' One row of the shared fit ledger: the polygon's geometry, the columns in
#' `fireSenseUtils::spreadFitAdditionalColNamesTxt`, `polygonID`, and, when the fit has an intercept, the
#' covariate centres (`fireSenseUtils::spreadFitCovCentreTxt`) and, when the fit simulated its members, the in-sample
#' fire sizes (`fireSenseUtils::spreadFitInSampleSizesTxt`) and, when the fit has fuel covariates, their median non-zero values
#' and sds (`fireSenseUtils::spreadFitCovMedianNonZeroTxt`, `fireSenseUtils::spreadFitCovSDTxt`). The `run` event writes it to the
#' ledger, and a held-out fold saves it beside its held-out result, so `fireSense_spreadPredict` reads
#' either the same way.
#'
#' @param sim a `simList`; uses `fireSenseUtils::escapeDefinitionOf()` (parameters `escapeSizeHa` and `fireCauses`), `.ELFind`, `studyArea`, `sppEquiv`, `nonForestedLCCGroups`,
#'   `missingLCCgroup`, `covMinMax_spread`, `covCentre_spread`, `climateCurvePars`, `covMedianNonZero_spread`, `covSD_spread` and `climateRef_spread`, and the parameters
#'   `spreadLink`, `spreadClimate`, `powerLawGamma`, `powerLawShape`, `lowerAsymptote`, `richardsNu` and `powerLawPcMargin`, and `powerLawCalibration`.
#' @param numIterations integer, the number of DEoptim generations.
#' @param objFunVal numeric, the objective values of the members in `params`.
#' @param params `data.table`, one row per member, with `hillSlope1` and `inflectionPoint1` (see `addFixedParsToLedger()`).
#' @param inSample `NULL` (a fold, or a fit without `simulateMembers`: no column) or `inSampleSizes()` of the fit.
#' @return an `sf` object with one row.
spreadFitLedgerRow <- function(sim, numIterations, objFunVal, params, inSample = NULL) {
  ## covMinMax_spread: prediction rescales covariates with it, exactly as this fit did
  df <- data.frame(I(list(numIterations)),
                   I(list(objFunVal)),
                   I(list(params)),
                   I(list(sim$sppEquiv)),
                   I(list(sim$nonForestedLCCGroups)),
                   I(list(sim$missingLCCgroup)),
                   I(list(sim$covMinMax_spread))) |>
    setNames(fireSenseUtils::spreadFitAdditionalColNamesTxt)
  ## covCentre_spread: a fit with an intercept centred its covariates, and prediction centres them alike.
  ## A fit without one has no such column, so its row is what it always was.
  if (!is.null(sim$covCentre_spread))
    df[[fireSenseUtils::spreadFitCovCentreTxt]] <- I(list(sim$covCentre_spread))
  ## climateCurvePars: a fit with curve terms (`<var>_curve`) stores `k` and `s` of each variable, and prediction builds the terms from them.
  ## A fit without curve terms has no such column, so its row is what it always was.
  if (!is.null(sim$climateCurvePars))
    df[[fireSenseUtils::spreadFitClimateCurveParsTxt]] <- I(list(sim$climateCurvePars))
  ## the fuel covariates' median non-zero value and sd, on the scale of the centres: what fireSenseUtils::spreadFuelEffect() reads.
  ## A fit without a fuel covariate has no such columns, and a row written before them has neither.
  if (!is.null(sim$covMedianNonZero_spread))
    df[[fireSenseUtils::spreadFitCovMedianNonZeroTxt]] <- I(list(sim$covMedianNonZero_spread))
  if (!is.null(sim$covSD_spread))
    df[[fireSenseUtils::spreadFitCovSDTxt]] <- I(list(sim$covSD_spread))
  ## inSample: the fit's own fires, simulated and observed (not a held-out fold's), for a report to read later
  if (!is.null(inSample))
    df[[fireSenseUtils::spreadFitInSampleSizesTxt]] <- I(list(inSample))
  ## spreadLink: a fit with another link than the default records it, with what the link needs (climateRef_spread
  ## for climate outside the ceiling, powerLawGamma for the power laws, lowerAsymptote and richardsNu for the links that
  ## have them -- never powerLawGamma for "richards", which does not use it), and prediction applies the same link.
  ## A fit with the default link writes none of these columns.
  ## spreadClimate: where climate enters, so prediction tells the activation from the ceiling (both have a power-law
  ## link, climateRef_spread and the p50 below). Like the link, only a fit with another than the default ("predictor")
  ## writes it; a row without it is read as "ceiling" when it has climateRef_spread and a power-law link, else "predictor".
  if (!identical(P(sim)$spreadClimate, "predictor")) df[[fireSenseUtils::spreadClimateTxt]] <- P(sim)$spreadClimate
  ## the activation's share of dud years, one per set, to judge the fit on (fireSenseUtils::activationDudYearShare())
  if (isActivation(P(sim)))
    df[[fireSenseUtils::spreadFitDudYearShareTxt]] <- I(list(fireSenseUtils::activationDudYearShare(params)))
  ## p_c fixed at p50: the row's maxAsymptote is p50, not a fitted value at its bound
  if (isTRUE(P(sim)$powerLawPcAtP50)) df[[fireSenseUtils::spreadFitPcAtP50Txt]] <- TRUE
  if (!identical(P(sim)$spreadLink, "ceiling")) df[[fireSenseUtils::spreadLinkTxt]] <- P(sim)$spreadLink
  if (!is.null(sim$climateRef_spread))
    df[[fireSenseUtils::spreadFitClimateRefTxt]] <- I(list(sim$climateRef_spread))
  if (P(sim)$spreadLink %in% fireSenseUtils::spreadAsymptoteLinks)
    df[[fireSenseUtils::spreadFitLowerAsymptoteTxt]] <- P(sim)$lowerAsymptote
  if (identical(P(sim)$spreadLink, "richards"))
    df[[fireSenseUtils::spreadFitRichardsNuTxt]] <- P(sim)$richardsNu
  if (P(sim)$spreadLink %in% fireSenseUtils::spreadPowerLawLinks) {
    df[[fireSenseUtils::spreadFitPowerLawGammaTxt]] <- powerLawGammaOf(sim)
    df[[fireSenseUtils::spreadFitPowerLawShapeTxt]] <- P(sim)$powerLawShape
    ## the calibration that set gamma (unless given) and p_c's upper bound (p50 + margin; p50 exactly with spreadClimate = "ceiling"), to compare p_c with
    df[[fireSenseUtils::spreadFitPowerLawPc50Txt]] <- sim$powerLawCalibration$p50
    df[[fireSenseUtils::spreadFitPowerLawR2Txt]] <- sim$powerLawCalibration$R2
    df[[fireSenseUtils::spreadFitPowerLawPcMarginTxt]] <- P(sim)$powerLawPcMargin
  }
  ## the escape definition the fit was made under: prediction and scoring check theirs against it
  ## (fireSenseUtils::checkEscapeDefinition())
  df <- data.frame(df, as.data.frame(fireSenseUtils::escapeDefinitionOf(sim)))
  ## the fuel scale of the covariates (fireSenseUtils::fuelLogToLinear() in spreadFitPrep): with the escape definition, what
  ## prediction must apply as this fit did, so fireSenseUtils::latestSpreadFits() keeps only rows that match the run
  df[[fireSenseUtils::spreadFitFuelScaleTxt]] <- fireSenseUtils::spreadFuelScale
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
