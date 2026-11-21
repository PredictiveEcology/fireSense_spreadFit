## spreadClimate = "activation": fuel sets each pixel's ceiling below p50 and climate is the linear predictor
## (fireSenseUtils::spreadLinkInputs(activationIntercept =)). The fit fits `activationIntercept` (after the covariates,
## before the sds) and not `maxAsymptote`; the ledger row has p50 in its place and records the mode.
## Helpers: helper-spreadLink.R.

pc50Col <- fireSenseUtils::spreadFitPowerLawPc50Txt
modeCol <- fireSenseUtils::spreadClimateTxt
actTxt <- fireSenseUtils::activationInterceptTxt
actParams <- list(spreadLink = "powerLawSmooth", spreadClimate = "activation")

test_that("'activation' fits activationIntercept right before yearSpreadSD and no maxAsymptote", {
  out <- fittedCO(actParams)
  a <- out$rec$deArgs
  p50 <- out$sim$powerLawCalibration$p50
  for (b in list(a$lower, a$upper)) {
    expect_false("maxAsymptote" %in% names(b))
    expect_identical(utils::tail(names(b), 2L), c(actTxt, fireSenseUtils::yearSpreadSDTxt))
  }
  expect_identical(unname(c(a$lower[[actTxt]], a$upper[[actTxt]])), fireSenseUtils::spreadActivationInterceptBounds)
  ## the objective learns it from the bounds, as it does the random effects
  expect_true(fireSenseUtils::spreadRandomEffectArgs(a$lower)$fitActivationIntercept)
  expect_named(out$sim$climateRef_spread, "CMDsm")
  expect_identical(a$climateRef, out$sim$climateRef_spread)
  expect_identical(a$powerLawPc50, p50)
  expect_identical(unname(a$lower[["CMDsm"]]), -50)   # the climate sign is free by default (coefSigns)
  ## its own bounds
  b2 <- fittedCO(c(actParams, list(activationInterceptBounds = c(-3, 4))))$rec$deArgs
  expect_identical(unname(c(b2$lower[[actTxt]], b2$upper[[actTxt]])), c(-3, 4))
})

test_that("'activation' writes the mode, p50 as maxAsymptote and the intercept to the ledger row", {
  out <- fittedCO(actParams)
  row <- out$rec$geoArgs$studyAreaFireSense
  p50 <- out$sim$powerLawCalibration$p50
  expect_identical(row[[modeCol]], "activation")
  expect_identical(row[[pc50Col]], p50)
  expect_identical(row[[refCol]][[1]], out$sim$climateRef_spread)
  pars <- row$params[[1]]
  expect_identical(names(pars)[1:3], c("maxAsymptote", "hillSlope1", "inflectionPoint1"))
  expect_true(all(pars$maxAsymptote == p50))
  expect_identical(utils::tail(names(pars), 2L), c(actTxt, fireSenseUtils::yearSpreadSDTxt))
})

test_that("'predictor' and 'ceiling' fit as before; only a ceiling row records its mode", {
  for (mode in c("predictor", "ceiling")) {
    out <- fittedCO(list(spreadLink = "powerLawSmooth", spreadClimate = mode))
    a <- out$rec$deArgs
    expect_identical(names(a$lower)[1], "maxAsymptote")
    expect_false(actTxt %in% names(a$lower))
    expect_null(fireSenseUtils::spreadRandomEffectArgs(a$lower)$fitActivationIntercept)
    row <- out$rec$geoArgs$studyAreaFireSense
    expect_false(actTxt %in% names(row$params[[1]]))
    ## the default writes no column, as before; fireSense_spreadPredict reads a row without it as "predictor"
    if (mode == "predictor") expect_false(modeCol %in% names(row)) else expect_identical(row[[modeCol]], mode)
  }
  expect_false(modeCol %in% names(fittedCO()$rec$geoArgs$studyAreaFireSense))
})

test_that("'activation' needs powerLawSmooth at the floor and a climate covariate", {
  for (p in list(list(spreadLink = "powerLaw"), list(lowerAsymptote = 0.05))) {
    s <- toySim(utils::modifyList(actParams, p), toyObjectsCO())
    mockFitAndLedger(s)
    expect_error(suppressMessages(SpaDES.core::spades(s)), "lowerAsymptote = lowerSpreadProb")
  }
  s <- toySim(actParams, toyObjectsCO("~ 0 + youngAge + class1 + class2 + nf"))
  mockFitAndLedger(s)
  expect_error(suppressMessages(SpaDES.core::spades(s)), "spreadClimate = 'activation' needs a climate covariate")
})

test_that("addFixedParsToLedger() puts maxAsymptote first only when given", {
  pb <- data.table::data.table(CMDsm = 1, nf = 2, activationIntercept = 0.5, yearSpreadSD = 0.3)
  got <- addFixedParsToLedger(pb, 0.27)
  expect_identical(names(got), c("maxAsymptote", "hillSlope1", "inflectionPoint1", names(pb)))
  expect_identical(got$maxAsymptote, 0.27)
  old <- data.table::data.table(maxAsymptote = 0.2, CMDsm = 1)
  expect_identical(addFixedParsToLedger(old), addFixedParsToLedger(old, NULL))
  expect_identical(names(addFixedParsToLedger(old)), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "CMDsm"))
})

test_that("the start-population gates, the debug call and the validation figures get p50", {
  out <- fittedCO(actParams)
  expect_identical(spreadLinkArgs(out$sim)$powerLawPc50, out$sim$powerLawCalibration$p50)
  fig <- fitWithPlots(c(actParams, list(.plots = "png")), toyObjectsCO())
  args <- fig$rec$data[[1]]$args
  expect_identical(args$powerLawPc50, fig$sim$powerLawCalibration$p50)
  expect_identical(args$climateRef, fig$sim$climateRef_spread)
})

test_that("with an intercept, 'activation' bounds it by activationFuelInterceptBounds; the other modes keep theirs", {
  formI <- "~ 1 + CMDsm + youngAge + class1 + class2 + nf"
  ic <- fireSenseUtils::spreadInterceptTxt
  a <- fittedCO(actParams, formula = formI)$rec$deArgs
  expect_identical(unname(c(a$lower[[ic]], a$upper[[ic]])), c(-6, 6))
  b <- fittedCO(c(actParams, list(activationFuelInterceptBounds = c(-3, 2))), formula = formI)$rec$deArgs
  expect_identical(unname(c(b$lower[[ic]], b$upper[[ic]])), c(-3, 2))
  for (mode in c("predictor", "ceiling")) {
    o <- fittedCO(list(spreadLink = "powerLawSmooth", spreadClimate = mode), formula = formI)
    v <- SpaDES.core::params(o$sim)$fireSense_spreadFit$upperAndLowerVal
    expect_identical(unname(c(o$rec$deArgs$lower[[ic]], o$rec$deArgs$upper[[ic]])), c(-v, v))
  }
})

test_that("an activation row stores each set's share of dud years; other rows have no such column", {
  out <- fittedCO(actParams)
  row <- out$rec$geoArgs$studyAreaFireSense
  share <- row[[fireSenseUtils::spreadFitDudYearShareTxt]][[1]]
  expect_identical(share, fireSenseUtils::activationDudYearShare(row$params[[1]]))
  expect_length(share, NROW(row$params[[1]]))
  expect_true(all(share >= 0 & share <= 1))
  expect_false(fireSenseUtils::spreadFitDudYearShareTxt %in%
                 names(fittedCO(list(spreadLink = "powerLawSmooth"))$rec$geoArgs$studyAreaFireSense))
})
