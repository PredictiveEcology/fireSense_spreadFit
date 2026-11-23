## The fuel covariates' median non-zero value and sd are stored in the ledger row, next to the centres, so the change in spread
## probability a fuel makes (fireSenseUtils::spreadFuelEffect()) can be found from the row alone. A fit without a fuel
## covariate, and a row written before the columns, have none.

fuelForm <- "~ 1 + CMDsm + youngAge + dom_agb_Pice_mar + dom_agb_Popu_tre + nf"
## the toy fuel columns named as the real ones are
toyObjectsFuel <- function(formula = fuelForm) {
  o <- toyObjectsI(formula)
  nonAnn <- o$fireSense_nonAnnualSpreadFitCovariates[[1]]
  data.table::setnames(nonAnn, c("class1", "class2"), c("dom_agb_Pice_mar", "dom_agb_Popu_tre"))
  o$sppEquiv <- data.table::data.table(LandR = c("Pice_mar", "Popu_tre"), FuelClass = c("Pice_mar", "Popu_tre"))
  o
}
fittedFuel <- function(formula = fuelForm) {
  rec <- new.env()
  sim <- toySim(list(), toyObjectsFuel(formula))
  mockFitAndLedger(sim, rec)
  list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
}

test_that("a fit with fuel covariates stores their median non-zero value and sd, as the data have them", {
  out <- fittedFuel()
  med <- out$sim$covMedianNonZero_spread
  expect_named(med, c("dom_agb_Pice_mar", "dom_agb_Popu_tre"))
  ## by hand: biomass over 1e4 (fireSenseUtils::fuelLinearRange), the toy years' pixels, 0 where the class is absent
  pice <- c(0, 300, 1000, 5000, 7803, 22601, 300, 1000, 5000) / 1e4
  popu <- c(40, 0, 2000, 0, 12000, 0, 500, 0, 9000) / 1e4
  ## the fit clears fuel on young pixels (mutually exclusive with youngAge), so those are not fuel in the data
  young <- c(0, 1, 0, 0, 0, 1, 0, 0, 1) == 1
  fuelInData <- function(x) { x[young] <- 0; x[x > 0] }
  expect_equal(med$dom_agb_Pice_mar, median(fuelInData(pice)), tolerance = 1e-3)
  expect_equal(med$dom_agb_Popu_tre, median(fuelInData(popu)), tolerance = 1e-3)
  expect_named(out$sim$covSD_spread, names(med))
  expect_true(all(unlist(out$sim$covSD_spread) > 0))
  ## the ledger row has both, beside the centres
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_identical(row[[fireSenseUtils::spreadFitCovMedianNonZeroTxt]][[1]], med)
  expect_identical(row[[fireSenseUtils::spreadFitCovSDTxt]][[1]], out$sim$covSD_spread)
  expect_true(fireSenseUtils::spreadFitCovCentreTxt %in% names(row))
  ## and fireSenseUtils reads the row: a change in probability per fuel, in 0-1
  effect <- fireSenseUtils::spreadFuelEffect(as.data.frame(row)[1, ])
  expect_named(effect, names(med))
  expect_false(anyNA(effect))
})

test_that("a fit without a fuel covariate adds no fuel columns to the ledger row", {
  out <- fittedI(formIcpt0 <- "~ 1 + CMDsm + youngAge + nf")
  expect_null(out$sim$covMedianNonZero_spread)
  expect_null(out$sim$covSD_spread)
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_false(any(c(fireSenseUtils::spreadFitCovMedianNonZeroTxt, fireSenseUtils::spreadFitCovSDTxt) %in% names(row)))
})

test_that("the new outputs are declared, so a cached event restores them", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_true(all(c("covMedianNonZero_spread", "covSD_spread") %in% md$outputObjects$objectName))
})
