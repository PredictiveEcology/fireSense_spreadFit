## spreadLink = "powerLaw": p = p_c - (p_c - lowerSpreadProb) * exp(-x / powerLawGamma) (fireSenseUtils::powerLawLink()).
## p_c is the fitted maxAsymptote, bounded by fireSenseUtils::spreadPowerLawPcBounds; gamma is a parameter, fixed. The
## link goes to the fit, the threshold calibration, the held-out simulation and the ledger row, which records it with
## gamma. With the default "ceiling" nothing changes. Helpers: helper-spreadLink.R.

gammaCol <- fireSenseUtils::spreadFitPowerLawGammaTxt

test_that("the defaults leave the fit, the ledger row and the cache keys as they were", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(md$parameters$default[[which(md$parameters$paramName == "powerLawGamma")]], 1.9)
  expect_identical(md$parameters$default[[which(md$parameters$paramName == "powerLawGamma")]],
                   fireSenseUtils::spreadPowerLawGamma)
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

test_that("with powerLaw the fit gets the link and gamma, p_c is bounded above the floor, other bounds are unchanged", {
  out <- fittedCO(list(spreadLink = "powerLaw", powerLawGamma = 2.5))
  a <- out$rec$deArgs
  expect_identical(a$link, "powerLaw")
  expect_identical(a$powerLawGamma, 2.5)
  expect_null(a$climateRef)
  pc <- fireSenseUtils::spreadPowerLawPcBounds
  expect_identical(c(unname(a$lower[["maxAsymptote"]]), unname(a$upper[["maxAsymptote"]])), pc)
  expect_gt(pc[1], fireSenseUtils::spreadProbFloor)
  ## the covariates' bounds are the ones of the default link: wide, sign only
  base <- fittedCO()$rec$deArgs
  expect_identical(a$lower[-1], base$lower[-1])
  expect_identical(a$upper[-1], base$upper[-1])
  expect_identical(unname(a$lower[["CMDsm"]]), 0)
  expect_identical(unname(a$upper[["youngAge"]]), 0)
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
})

test_that("the held-out simulation and the fold's ledger row use the fit's link and gamma", {
  run <- function(params) {
    rec <- new.env(); rec$fnArgs <- list()
    sim <- toySim(c(list(heldOutFold = 1L, simulateMembers = 2L, doFit = FALSE), params), toyObjectsCO())
    mockFitAndLedger(sim, rec)
    mockInModule(sim, Cache = function(FUN, pop, fnArgs, ...) {
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
  expect_identical(maxAsymptoteBoundsOf("powerLaw"), fireSenseUtils::spreadPowerLawPcBounds)
  pl <- maxAsymptoteBoundsOf("powerLaw")
  expect_identical(estimateSpreadParams(form, annual, "lower", 9, maxAsymptoteBounds = pl),
                   c(maxAsymptote = pl[1], CMDsm = 0, youngAge = -9))
  expect_identical(estimateSpreadParams(form, annual, "upper", 9, maxAsymptoteBounds = pl),
                   c(maxAsymptote = pl[2], CMDsm = 9, youngAge = 0))
})
