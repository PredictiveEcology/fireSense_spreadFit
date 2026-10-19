#' Calibrate the `"powerLaw"` link on the ELF's flammable map
#'
#' With `spreadLink = "powerLaw"`, runs `fireSenseUtils::calibratePowerLaw()` once per ELF, on `flammableRTM`
#' (the latest vegetation year's flammable map) and with the fit's own escape size and jump settings, and
#' stores the result as `sim$powerLawCalibration`. It always runs, because its `p50` bounds `p_c`; a numeric
#' `powerLawGamma` replaces only `gamma` ([powerLawGammaOf()]). The `Cache()` key is the map and the arguments
#' (not `cores`: the result does not depend on it) and the bodies of the calibration and its callee.
#'
#' @param sim a `simList`.
#' @return the `simList`; unchanged unless `spreadLink` is `"powerLaw"`.
calibratePowerLawFor <- function(sim) {
  if (!identical(P(sim)$spreadLink, "powerLaw")) return(sim)
  if (!inherits(sim$flammableRTM, "SpatRaster"))
    stop("fireSense_spreadFit: spreadLink = 'powerLaw' calibrates on `flammableRTM`, a SpatRaster (from fireSense_dataPrepFit); ",
         "supply it, or use spreadLink = 'ceiling'")
  sim$powerLawCalibration <- Cache(
    fireSenseUtils::calibratePowerLaw,
    flammable = sim$flammableRTM, pGrid = P(sim)$powerLawPGrid, nFires = P(sim)$powerLawFires,
    cores = calibrationCores(P(sim)$cores), seed = .elfSeed(sim$.ELFind),
    escapeSizeHa = P(sim)$escapeSizeHa, jumpTries = P(sim)$jumpTries, jumpMeanDist = P(sim)$jumpMeanDist,
    omitArgs = "cores",
    .cacheExtra = list(calibratePowerLaw = deparse(fireSenseUtils::calibratePowerLaw),
                       escapeSizePixels = deparse(fireSenseUtils::escapeSizePixels)),
    useCache = P(sim)$useCache_DE)
  sim
}

## Local forked processes for the calibration: the module's `cores` if it is one number, else (host names) as
## many as it lists; never more than this machine has
calibrationCores <- function(cores) {
  n <- if (is.numeric(cores) && length(cores) == 1L && !is.na(cores)) cores else length(cores)
  max(1L, min(as.integer(n), parallel::detectCores(), na.rm = TRUE))
}

#' The `"powerLaw"` link's `gamma` for this fit
#'
#' `powerLawGamma` if it is a number, else the calibrated one; `fireSenseUtils::spreadPowerLawGamma` when there is
#' no calibration (the other links, which do not use it).
#'
#' @param sim a `simList`.
#' @return A number.
powerLawGammaOf <- function(sim) {
  g <- P(sim)$powerLawGamma
  if (!is.na(g)) g else if (!is.null(sim$powerLawCalibration)) sim$powerLawCalibration$gamma else fireSenseUtils::spreadPowerLawGamma
}
