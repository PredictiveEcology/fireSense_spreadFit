## spreadLink = "powerLaw": p = p_c - (p_c - lowerSpreadProb) * exp(-x / powerLawGamma) (fireSenseUtils::powerLawLink()).
## p_c is the fitted maxAsymptote, bounded above per ELF by the calibrated p_c50 + powerLawPcMargin; gamma is calibrated on
## the ELF's flammable map unless powerLawGamma is a number (fixed, not fitted). The
## link goes to the fit, the threshold calibration, the held-out simulation and the ledger row, which records it with
## gamma. With the default "ceiling" nothing changes. Helpers: helper-spreadLink.R.

gammaCol <- fireSenseUtils::spreadFitPowerLawGammaTxt

test_that("the defaults leave the fit, the ledger row and the cache keys as they were", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_true("powerLaw" %in% fireSenseUtils::spreadLinkOptions)
  out <- fittedCO()
  expect_null(out$rec$deArgs$link)
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_false(any(c(linkCol, gammaCol) %in% names(row)))
  expect_identical(unname(out$rec$deArgs$upper[["maxAsymptote"]]), fireSenseUtils::spreadProbCeiling)
  expect_identical(unname(out$rec$deArgs$lower[["maxAsymptote"]]), fireSenseUtils::spreadProbFloor)
  off <- omitArgsCO()
  expect_true(all(vapply(off$omit, function(o) "powerLawGamma" %in% o, logical(1))))
  expect_null(off$rswd$link)
})

test_that("the defaults are the powerLaw link, gamma calibrated, no runaway penalty", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  dflt <- function(p) md$parameters$default[[which(md$parameters$paramName == p)]]
  expect_identical(dflt("spreadLink"), "powerLaw")
  expect_true(is.na(dflt("powerLawGamma")))
  expect_identical(dflt("powerLawPcMargin"), 0.005)
  expect_false(dflt("penaliseRunaways"))
  ## a run with the module's own defaults (only the calibration made small) fits the powerLaw link, calibrated
  out <- fittedCO(list(spreadLink = NULL, penaliseRunaways = NULL))
  expect_identical(out$rec$deArgs$link, "powerLaw")
  expect_false(out$rec$deArgs$penaliseRunaways)
  expect_identical(out$rec$geoArgs$studyAreaFireSense[[linkCol]], "powerLaw")
  expect_identical(out$rec$deArgs$powerLawGamma, out$sim$powerLawCalibration$gamma)
})

test_that("powerLaw calibrates gamma on flammableRTM, bounds p_c by the calibrated p_c50, and the ledger row records them", {
  out <- fittedCO(list(spreadLink = "powerLaw"))
  a <- out$rec$deArgs
  cal <- out$sim$powerLawCalibration
  expect_identical(a$link, "powerLaw")
  expect_identical(a$powerLawGamma, cal$gamma)
  expect_gt(cal$gamma, 0)
  ## the calibration is the package function on the same map and settings
  again <- fireSenseUtils::calibratePowerLaw(calFlammable(), pGrid = calParams$powerLawPGrid, nFires = 20L,
                                             seed = .elfSeed("9.9"),
                                             escapeSizeHa = 50, jumpTries = 20, jumpMeanDist = 3)
  expect_identical(cal$gamma, again$gamma)
  ## DEoptim's bounds for p_c come from it
  expect_equal(c(unname(a$lower[["maxAsymptote"]]), unname(a$upper[["maxAsymptote"]])),
               c(fireSenseUtils::spreadProbFloor + 0.01, cal$p50 + 0.005))
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_identical(row[[gammaCol]], cal$gamma)
  expect_identical(row[[fireSenseUtils::spreadFitPowerLawPc50Txt]], cal$p50)
  expect_identical(row[[fireSenseUtils::spreadFitPowerLawR2Txt]], cal$R2)
  expect_identical(row[[fireSenseUtils::spreadFitPowerLawPcMarginTxt]], 0.005)
  ## the margin is a parameter
  wide <- fittedCO(list(spreadLink = "powerLaw", powerLawPcMargin = 0.01))
  expect_equal(unname(wide$rec$deArgs$upper[["maxAsymptote"]]), wide$sim$powerLawCalibration$p50 + 0.01)
})

test_that("a numeric powerLawGamma overrides gamma only: the calibration still runs and bounds p_c", {
  out <- fittedCO(list(spreadLink = "powerLaw", powerLawGamma = 2.5))
  a <- out$rec$deArgs
  expect_identical(a$link, "powerLaw")
  expect_identical(a$powerLawGamma, 2.5)
  expect_false(identical(out$sim$powerLawCalibration$gamma, 2.5))
  expect_equal(unname(a$upper[["maxAsymptote"]]), out$sim$powerLawCalibration$p50 + 0.005)
  expect_identical(out$rec$geoArgs$studyAreaFireSense[[gammaCol]], 2.5)
  expect_null(a$climateRef)
  expect_gt(a$lower[["maxAsymptote"]], fireSenseUtils::spreadProbFloor)
  ## the covariates' bounds are the ones of the default link: wide, sign only
  base <- fittedCO()$rec$deArgs
  expect_identical(a$lower[-1], base$lower[-1])
  expect_identical(a$upper[-1], base$upper[-1])
  expect_identical(unname(a$lower[["CMDsm"]]), 0)
  expect_identical(unname(a$upper[["youngAge"]]), 0)
})

test_that("ceiling and climateOutside run no calibration and need no flammableRTM", {
  for (lk in c("ceiling", "climateOutside")) {
    rec <- new.env()
    sim <- toySim(list(spreadLink = lk), c(toyObjectsCO(), list(flammableRTM = NULL)))
    mockFitAndLedger(sim, rec)
    sim <- suppressMessages(SpaDES.core::spades(sim))
    expect_null(sim$powerLawCalibration)
  }
  s <- toySim(list(spreadLink = "powerLaw"), c(toyObjectsCO(), list(flammableRTM = NULL)))
  mockFitAndLedger(s)
  expect_error(suppressMessages(SpaDES.core::spades(s)), "flammableRTM")
})

test_that("the ledger row records the link and gamma, and no climate reference", {
  out <- fittedCO(list(spreadLink = "powerLaw", powerLawGamma = 2.5))
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_identical(row[[linkCol]], "powerLaw")
  expect_identical(row[[gammaCol]], 2.5)
  expect_false(refCol %in% names(row))
  expect_identical(names(row$params[[1]])[1:4], c("maxAsymptote", "hillSlope1", "inflectionPoint1", "CMDsm"))
  ## climateOutside is unchanged: its link, its reference, no gamma
  co <- fittedCO(list(spreadLink = "climateOutside"))$rec$geoArgs$studyAreaFireSense
  expect_identical(co[[linkCol]], "climateOutside")
  expect_true(refCol %in% names(co))
  expect_false(gammaCol %in% names(co))
})

test_that("an unknown spreadLink, powerLaw with the upper-tail link, and a bad gamma stop", {
  run <- function(params) { sim <- toySim(params, toyObjectsCO()); mockFitAndLedger(sim); suppressMessages(SpaDES.core::spades(sim)) }
  expect_error(run(list(spreadLink = "power")), "spreadLink must be one of")
  expect_error(run(list(spreadLink = "powerLaw", link = "logistic3pUpper")), "cannot be combined")
  expect_error(run(list(spreadLink = "powerLaw", powerLawGamma = 0)), "powerLawGamma must be one positive number")
  expect_error(run(list(spreadLink = "powerLaw", powerLawGamma = -1)), "powerLawGamma must be one positive number")
  expect_error(run(list(spreadLink = "powerLaw", powerLawGamma = c(1, 2))), "powerLawGamma must be one positive number")
})

test_that("the threshold calibration is given the link and gamma, and gamma is in the cache keys only with powerLaw", {
  on <- omitArgsCO(list(spreadLink = "powerLaw", powerLawGamma = 2.5))
  expect_false(any(vapply(on$omit, function(o) "powerLawGamma" %in% o, logical(1))))
  expect_identical(on$rswd$link, "powerLaw")
  expect_identical(on$rswd$powerLawGamma, 2.5)
  src <- paste(deparse(body(runSpreadWithoutDEoptim)), collapse = "\n")
  ## the debug objective call, the gates, and the trials
  expect_equal(lengths(regmatches(src, gregexpr("powerLawGamma = powerLawGamma", src))), 3L)
  expect_identical(formals(runSpreadWithoutDEoptim)$powerLawGamma, quote(fireSenseUtils::spreadPowerLawGamma))
  ## with gamma calibrated, the calibrated one goes to the threshold calibration
  cal <- omitArgsCO(list(spreadLink = "powerLaw"))
  expect_gt(cal$rswd$powerLawGamma, 0)
  expect_false(isTRUE(all.equal(cal$rswd$powerLawGamma, 2.5)))
})

test_that("the held-out simulation and the fold's ledger row use the fit's link and gamma", {
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
  on <- run(list(spreadLink = "powerLaw", powerLawGamma = 2.5))
  expect_identical(on$rec$fnArgs$link, "powerLaw")
  expect_identical(on$rec$fnArgs$powerLawGamma, 2.5)
  expect_identical(on$sim$spreadFitHeldOut$fit[[linkCol]], "powerLaw")
  expect_identical(on$sim$spreadFitHeldOut$fit[[gammaCol]], 2.5)
  off <- run(list())
  expect_null(off$rec$fnArgs$link)
  expect_false("powerLawGamma" %in% names(off$rec$fnArgs))
  expect_false(gammaCol %in% names(off$sim$spreadFitHeldOut$fit))
})

test_that("maxAsymptoteBoundsOf() and estimateSpreadParams(): p_c's bounds follow the link, nothing else moves", {
  form <- "~ 0 + CMDsm + youngAge"
  annual <- list(year2001 = data.table::data.table(pixelID = 1L, CMDsm = 1, youngAge = 0))
  expect_identical(maxAsymptoteBoundsOf(), c(fireSenseUtils::spreadProbFloor, fireSenseUtils::spreadProbCeiling))
  expect_identical(maxAsymptoteBoundsOf("climateOutside"), maxAsymptoteBoundsOf("ceiling"))
  pl <- maxAsymptoteBoundsOf("powerLaw", p50 = 0.26, margin = 0.005)
  expect_equal(pl, c(fireSenseUtils::spreadProbFloor + 0.01, 0.265))
  expect_error(maxAsymptoteBoundsOf("powerLaw", p50 = 0.5), "outside")
  expect_identical(estimateSpreadParams(form, annual, "lower", 9, maxAsymptoteBounds = pl),
                   c(maxAsymptote = pl[1], CMDsm = 0, youngAge = -9))
  expect_identical(estimateSpreadParams(form, annual, "upper", 9, maxAsymptoteBounds = pl),
                   c(maxAsymptote = pl[2], CMDsm = 9, youngAge = 0))
})

## calibratePowerLaw() forks one process per fire; only a job that fits (or validates) uses its result
calibrationCalls <- function(params, objects = list()) {
  rec <- new.env()
  rec$n <- 0L
  sim <- toySim(c(list(spreadLink = "powerLaw"), params), c(toyObjectsCO(), objects))
  mockFitAndLedger(sim, rec)
  mockInModule(sim, Cache = function(FUN, ...) {
    if (!identical(FUN, fireSenseUtils::calibratePowerLaw)) return(reproducible::Cache(FUN, ...))
    rec$n <- rec$n + 1L
    runCalibration(FUN, ...)
  })
  sim <- suppressMessages(SpaDES.core::spades(sim))
  rec$n
}

test_that("calibratePowerLaw runs when the job fits, and not when it uses the ledger's row", {
  row <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  expect_identical(calibrationCalls(list()), 1L)                              # no row: fit
  expect_identical(calibrationCalls(list(doFit = TRUE), row), 1L)             # refit over a row
  expect_identical(calibrationCalls(list(doFit = FALSE), row), 0L)            # hindcast: use the row
  expect_identical(calibrationCalls(list(doFit = NA), row), 0L)               # row found, reused
})
