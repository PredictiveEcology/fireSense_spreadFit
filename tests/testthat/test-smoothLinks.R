## spreadLink = "powerLawSmooth" and "richards": links with a lower asymptote `lowerAsymptote` (0 to spreadProbFloor), from
## fireSenseUtils. "powerLawSmooth" is calibrated like "powerLaw" (same gamma and bounds of p_c); "richards" is never
## calibrated (no flammableRTM needed, no powerLawGamma used or written). The ledger row records the link and its
## parameters. Helpers: helper-spreadLink.R, helper-toyInputs.R.

lowerCol <- fireSenseUtils::spreadFitLowerAsymptoteTxt
nuCol <- fireSenseUtils::spreadFitRichardsNuTxt
gammaCol <- fireSenseUtils::spreadFitPowerLawGammaTxt
shapeCol <- fireSenseUtils::spreadFitPowerLawShapeTxt

test_that("the parameters default to the floor and 0.5, and the new links are options", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  dflt <- function(p) md$parameters$default[[which(md$parameters$paramName == p)]]
  expect_identical(dflt("lowerAsymptote"), fireSenseUtils::spreadProbFloor)
  expect_identical(dflt("richardsNu"), 0.5)
  expect_true(all(c("powerLawSmooth", "richards") %in% fireSenseUtils::spreadLinkOptions))
  expect_identical(lowerCol, "lowerAsymptote_spread")
  expect_identical(nuCol, "richardsNu_spread")
})

test_that("powerLawSmooth is calibrated like powerLaw; the link, lowerAsymptote and the objective's arguments reach the fit and the ledger row", {
  out <- fittedCO(list(spreadLink = "powerLawSmooth", lowerAsymptote = 0.05))
  a <- out$rec$deArgs
  cal <- out$sim$powerLawCalibration
  expect_identical(a$link, "powerLawSmooth")
  expect_identical(a$lowerAsymptote, 0.05)
  expect_false("richardsNu" %in% names(a))
  expect_identical(a$powerLawGamma, cal$gamma)
  ## the bounds of p_c are the powerLaw's
  expect_equal(c(unname(a$lower[["maxAsymptote"]]), unname(a$upper[["maxAsymptote"]])),
               c(fireSenseUtils::spreadProbFloor + 0.01, cal$p50 + 0.005))
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_identical(row[[linkCol]], "powerLawSmooth")
  expect_identical(row[[lowerCol]], 0.05)
  expect_identical(row[[gammaCol]], cal$gamma)
  expect_identical(row[[shapeCol]], 1)
  expect_false(nuCol %in% names(row))
  ## its key has the lower asymptote and not nu
  om <- omitArgsCO(list(spreadLink = "powerLawSmooth"))$omit[[1]]
  expect_false("lowerAsymptote" %in% om)
  expect_true("richardsNu" %in% om)
})

test_that("richards runs no calibration, needs no flammableRTM, and writes no gamma", {
  rec <- new.env()
  sim <- toySim(list(spreadLink = "richards", lowerAsymptote = 0, richardsNu = 2, powerLawGamma = 5),
                c(toyObjectsCO(), list(flammableRTM = NULL)))
  mockFitAndLedger(sim, rec)
  calls <- 0L
  mockInModule(sim, Cache = function(FUN, ...) {
    if (identical(FUN, fireSenseUtils::calibratePowerLaw)) calls <<- calls + 1L
    reproducible::Cache(FUN, ...)
  })
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_identical(calls, 0L)
  expect_null(sim$powerLawCalibration)
  a <- rec$deArgs
  expect_identical(a$link, "richards")
  expect_identical(a$lowerAsymptote, 0)
  expect_identical(a$richardsNu, 2)
  expect_equal(c(unname(a$lower[["maxAsymptote"]]), unname(a$upper[["maxAsymptote"]])), fireSenseUtils::spreadPowerLawPcBounds)
  row <- rec$geoArgs$studyAreaFireSense
  expect_identical(row[[linkCol]], "richards")
  expect_identical(row[[lowerCol]], 0)
  expect_identical(row[[nuCol]], 2)
  expect_false(any(c(gammaCol, shapeCol, fireSenseUtils::spreadFitPowerLawPc50Txt) %in% names(row)))
  ## the key has the lower asymptote and nu, and not gamma or shape
  om <- omitArgsCO(list(spreadLink = "richards"))$omit[[1]]
  expect_false(any(c("lowerAsymptote", "richardsNu") %in% om))
  expect_true(all(c("powerLawGamma", "powerLawShape") %in% om))
  ## the power-law links still need the map
  s <- toySim(list(spreadLink = "powerLawSmooth"), c(toyObjectsCO(), list(flammableRTM = NULL)))
  mockFitAndLedger(s)
  expect_error(suppressMessages(SpaDES.core::spades(s)), "flammableRTM")
})

test_that("powerLaw and the default link are as they were: no new ledger columns, the new arguments left out of the key", {
  for (lk in c("powerLaw", "ceiling")) {
    out <- fittedCO(list(spreadLink = lk))
    expect_false(any(c(lowerCol, nuCol) %in% names(out$rec$geoArgs$studyAreaFireSense)))
    expect_false(any(c("lowerAsymptote", "richardsNu") %in% names(out$rec$deArgs)))
    expect_true(all(c("lowerAsymptote", "richardsNu") %in% omitArgsCO(list(spreadLink = lk))$omit[[1]]))
  }
})

test_that("a lowerAsymptote outside 0 to the floor, or a negative richardsNu, stops", {
  run <- function(params) { sim <- toySim(params, toyObjectsCO()); mockFitAndLedger(sim); suppressMessages(SpaDES.core::spades(sim)) }
  for (lk in c("powerLawSmooth", "richards")) {
    expect_error(run(list(spreadLink = lk, lowerAsymptote = 0.2)), "lowerAsymptote must be")
    expect_error(run(list(spreadLink = lk, lowerAsymptote = -0.01)), "lowerAsymptote must be")
    expect_error(run(list(spreadLink = lk, lowerAsymptote = c(0, 0.1))), "lowerAsymptote must be")
  }
  expect_error(run(list(spreadLink = "richards", richardsNu = -1)), "richardsNu must be")
  expect_error(run(list(spreadLink = "richards", link = "logistic3pUpper")), "cannot be combined")
  ## ignored with the other links
  expect_no_error(run(list(spreadLink = "powerLaw", lowerAsymptote = 0.5, richardsNu = -1)))
})

test_that("the held-out simulation gets the link's arguments, and its ledger row the link's columns", {
  run <- function(params) {
    rec <- new.env(); rec$fnArgs <- list()
    sim <- toySim(c(list(heldOutFold = 1L, simulateMembers = 2L, doFit = FALSE), params), toyObjectsCO())
    mockFitAndLedger(sim, rec)
    mockInModule(sim, Cache = function(FUN, pop, fnArgs, ...) {
      if (identical(FUN, fireSenseUtils::calibratePowerLaw)) return(runCalibration(FUN, ...))
      if (!is.function(FUN)) return(FUN)
      rec$fnArgs <- fnArgs
      data.table::data.table(member = 1L, yr = "year2001", rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
    list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
  }
  r <- run(list(spreadLink = "richards", lowerAsymptote = 0.02, richardsNu = 3))
  expect_identical(r$rec$fnArgs$link, "richards")
  expect_identical(r$rec$fnArgs$lowerAsymptote, 0.02)
  expect_identical(r$rec$fnArgs$richardsNu, 3)
  expect_false("powerLawGamma" %in% names(r$rec$fnArgs))
  expect_identical(r$sim$spreadFitHeldOut$fit[[nuCol]], 3)
  expect_false(gammaCol %in% names(r$sim$spreadFitHeldOut$fit))
  s <- run(list(spreadLink = "powerLawSmooth", lowerAsymptote = 0.05, powerLawGamma = 2.5))
  expect_identical(s$rec$fnArgs$lowerAsymptote, 0.05)
  expect_identical(s$rec$fnArgs$powerLawGamma, 2.5)
  expect_false("richardsNu" %in% names(s$rec$fnArgs))
  off <- run(list(spreadLink = "powerLaw"))
  expect_false(any(c("lowerAsymptote", "richardsNu") %in% names(off$rec$fnArgs)))
})

test_that("the debug event and the wholeLandscapeSpread switch accept the new links", {
  for (lk in c("powerLawSmooth", "richards")) {
    d <- debugArgs(list(spreadLink = lk, lowerAsymptote = 0.07, richardsNu = 1))
    expect_identical(d$link, lk)
    expect_identical(d$lowerAsymptote, 0.07)
    expect_null(checkWholeLandscapeSpread(list(wholeLandscapeSpread = TRUE, spreadLink = lk, penaliseRunaways = FALSE)))
  }
  expect_error(checkWholeLandscapeSpread(list(wholeLandscapeSpread = TRUE, spreadLink = "ceiling")), "richards")
  expect_identical(maxAsymptoteBoundsOf("richards"), fireSenseUtils::spreadPowerLawPcBounds)
  expect_equal(maxAsymptoteBoundsOf("powerLawSmooth", p50 = 0.26, margin = 0.005), c(fireSenseUtils::spreadProbFloor + 0.01, 0.265))
})
