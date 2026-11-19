## spreadClimate = "ceiling": climate moves the power-law ceiling p_c toward the landscape's runaway threshold p50
## (fireSenseUtils::spreadLinkInputs()) instead of entering the linear predictor under a fixed ceiling. The fit gets
## the climate reference as climateOutside does, plus p50, and bounds p_c above by p50 exactly. With the default
## "predictor" nothing changes. Helpers: helper-spreadLink.R.

pc50Col <- fireSenseUtils::spreadFitPowerLawPc50Txt

test_that("'predictor' (the default) leaves the fit, the ledger row and the cache keys as they were", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(md$parameters$default[[which(md$parameters$paramName == "spreadClimate")]], "predictor")
  out <- fittedCO(list(spreadLink = "powerLaw"))
  expect_null(out$sim$climateRef_spread)
  expect_false("powerLawPc50" %in% names(out$rec$deArgs))
  expect_false(refCol %in% names(out$rec$geoArgs$studyAreaFireSense))
  expect_equal(unname(out$rec$deArgs$upper[["maxAsymptote"]]), out$sim$powerLawCalibration$p50 + 0.005)
})

test_that("'ceiling' gives the fit the climate reference and p50, and bounds p_c above by p50", {
  out <- fittedCO(list(spreadLink = "powerLaw", spreadClimate = "ceiling"))
  a <- out$rec$deArgs
  p50 <- out$sim$powerLawCalibration$p50
  expect_named(out$sim$climateRef_spread, "CMDsm")
  expect_equal(out$sim$climateRef_spread$CMDsm, cmdMean(out$sim))
  expect_identical(a$climateRef, out$sim$climateRef_spread)
  expect_identical(a$powerLawPc50, p50)
  expect_identical(unname(a$upper[["maxAsymptote"]]), p50)
  expect_equal(unname(a$lower[["maxAsymptote"]]), fireSenseUtils::spreadProbFloor + 0.01)
  expect_identical(unname(a$lower[["CMDsm"]]), 0)   # the climate sign constraint is unchanged
  ## its cache key carries both
  on <- omitArgsCO(list(spreadLink = "powerLaw", spreadClimate = "ceiling"))
  expect_false(any(c("climateRef", "powerLawPc50") %in% on$omit[[1]]))
})

test_that("'ceiling' writes climateRef_spread, the link and p50 to the ledger row", {
  out <- fittedCO(list(spreadLink = "powerLawSmooth", spreadClimate = "ceiling"))
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_identical(row[[refCol]][[1]], out$sim$climateRef_spread)
  expect_identical(row[[linkCol]], "powerLawSmooth")
  expect_identical(row[[pc50Col]], out$sim$powerLawCalibration$p50)
})

test_that("'ceiling' needs a power-law link and a climate covariate", {
  expect_error(suppressMessages(SpaDES.core::spades({
    s <- toySim(list(spreadLink = "ceiling", spreadClimate = "ceiling"), toyObjectsCO()); mockFitAndLedger(s); s })),
    "needs spreadLink one of")
  s <- toySim(list(spreadLink = "powerLaw", spreadClimate = "ceiling"),
              toyObjectsCO("~ 0 + youngAge + class1 + class2 + nf"))
  mockFitAndLedger(s)
  expect_error(suppressMessages(SpaDES.core::spades(s)), "spreadClimate = 'ceiling' needs a climate covariate")
  s <- toySim(list(spreadLink = "powerLaw", spreadClimate = "wet"), toyObjectsCO())
  expect_error(suppressMessages(SpaDES.core::spades(s)), "spreadClimate must be")
})

test_that("the validation simulation behind the response curves gets the reference and p50", {
  out <- fitWithPlots(list(spreadLink = "powerLaw", spreadClimate = "ceiling", .plots = "png"), toyObjectsCO())
  args <- out$rec$data[[1]]$args
  expect_identical(args$powerLawPc50, out$sim$powerLawCalibration$p50)
  expect_identical(args$climateRef, out$sim$climateRef_spread)
})

test_that("the start-population gates and the debug objective call get p50", {
  out <- fittedCO(list(spreadLink = "powerLaw", spreadClimate = "ceiling"))
  expect_true("powerLawPc50" %in% names(formals(spreadGates)))
  expect_true("powerLawPc50" %in% names(formals(runSpreadWithoutDEoptim)))
  expect_identical(spreadLinkArgs(out$sim)$powerLawPc50, out$sim$powerLawCalibration$p50)
  expect_false("powerLawPc50" %in% names(spreadLinkArgs(fittedCO(list(spreadLink = "powerLaw"))$sim)))
})
