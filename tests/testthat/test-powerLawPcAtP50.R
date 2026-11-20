## powerLawPcAtP50 = TRUE fixes p_c (maxAsymptote) at the calibrated p50. As with spreadClimate = "activation", it is left
## out of the bounds and the objective puts p50 in its place (fireSenseUtils::.objfunSpreadFit(pcAtP50 = TRUE)): no
## zero-width bound, which the emulator's search space cannot scale. Helpers: helper-spreadLink.R, helper-seededRun.R.

test_that("maxAsymptoteBoundsOf(pcAtP50 = TRUE) leaves maxAsymptote out, still checking p50; FALSE is as before", {
  expect_null(maxAsymptoteBoundsOf("powerLawSmooth", p50 = 0.26, pcAtP50 = TRUE))
  expect_null(maxAsymptoteBoundsOf("powerLaw", p50 = 0.26, pcAtP50 = TRUE))
  expect_identical(maxAsymptoteBoundsOf("powerLaw", p50 = 0.26, pcAtP50 = FALSE), maxAsymptoteBoundsOf("powerLaw", p50 = 0.26))
  expect_error(maxAsymptoteBoundsOf("powerLaw", p50 = 0.5, pcAtP50 = TRUE), "outside")
})

test_that("powerLawPcAtP50 is off by default and, when TRUE, the fit has no maxAsymptote and gets p50 and pcAtP50", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_false(md$parameters$default[[which(md$parameters$paramName == "powerLawPcAtP50")]])
  out <- fittedCO(list(spreadLink = "powerLawSmooth", powerLawShape = 1, powerLawPcAtP50 = TRUE))
  a <- out$rec$deArgs
  p50 <- out$sim$powerLawCalibration$p50
  expect_false("maxAsymptote" %in% c(names(a$lower), names(a$upper)))
  expect_true(all(a$upper > a$lower))                                   # no zero-width bound
  expect_identical(unname(a$upper[["CMDsm"]]) > 0, TRUE)                # the other bounds are unchanged
  expect_identical(a$powerLawPc50, p50)
  expect_true(a$pcAtP50)
  expect_null(a$climateRef)                                             # still the "predictor" link
  ## the ledger row carries p50 as p_c, first, as the objective used it
  pars <- out$rec$geoArgs$studyAreaFireSense$params[[1]]
  expect_identical(names(pars)[1:3], c("maxAsymptote", "hillSlope1", "inflectionPoint1"))
  expect_true(all(pars$maxAsymptote == p50))
  ## and says it was fixed, not fitted at a bound
  expect_true(out$rec$geoArgs$studyAreaFireSense[[fireSenseUtils::spreadFitPcAtP50Txt]])
  off <- fittedCO(list(spreadLink = "powerLawSmooth", powerLawShape = 1))
  expect_false(fireSenseUtils::spreadFitPcAtP50Txt %in% names(off$rec$geoArgs$studyAreaFireSense))
  expect_lt(unname(off$rec$deArgs$lower[["maxAsymptote"]]), p50)
  expect_false(any(c("pcAtP50", "powerLawPc50") %in% names(off$rec$deArgs)))
})

test_that("powerLawPcAtP50 = TRUE stops with the ceiling and with a link without p_c50", {
  s <- toySim(list(spreadLink = "powerLaw", spreadClimate = "ceiling", powerLawPcAtP50 = TRUE), toyObjectsCO())
  mockFitAndLedger(s)
  expect_error(suppressMessages(SpaDES.core::spades(s)), "powerLawPcAtP50 = TRUE cannot be combined with spreadClimate")
  s <- toySim(list(spreadLink = "richards", powerLawPcAtP50 = TRUE), toyObjectsCO())
  mockFitAndLedger(s)
  expect_error(suppressMessages(SpaDES.core::spades(s)), "powerLawPcAtP50 = TRUE needs spreadLink")
})

test_that("an emulator start works with powerLawPcAtP50: finite search space, the objective told, a population in bounds", {
  rec <- seededRun(params = list(initialPopMethod = "emulator", spreadLink = "powerLawSmooth", powerLawPcAtP50 = TRUE),
                   mocks = emulatorMock)
  lo <- rec$lower; up <- rec$upper
  expect_false("maxAsymptote" %in% names(lo))
  ## the objective the emulator evaluates fixes p_c at p50, as the fit's does
  expect_true(rec$emu$objArgs$pcAtP50)
  expect_false(is.null(rec$emu$objArgs$powerLawPc50))
  ## the unit cube of these bounds is finite; #141's equal bounds would give NaN there
  sp <- emulatorSpace(lo, up)
  u <- sp$toUnit(rbind(lo + 0.5 * (up - lo)))
  expect_false(anyNA(u))
  pinned <- c(maxAsymptote = 0.26, lo)
  expect_true(anyNA(emulatorSpace(pinned, c(maxAsymptote = 0.26, up))$toUnit(rbind(pinned))))
  ## and the real search, on a toy objective over these bounds, returns a population inside them
  mid <- lo + 0.4 * (up - lo)
  value <- function(P, seed) withr::with_seed(seed, 100 + rowSums(sweep(P[, names(lo), drop = FALSE], 2, mid)^2) +
                                                     stats::rnorm(nrow(P), 0, 0.1))
  ctl <- list(designPerParam = 5L, globalMaxEvals = 60L, maxEvals = 80L, minEvals = 40L, batchSize = 8L,
              nCandidates = 200L, nPerturb = 100L, gpThreads = 1L)
  res <- suppressMessages(emulateObjective(value, lo, up, NP = 4L, control = ctl, seed = 1L))
  expect_identical(colnames(res$pop), names(lo))
  expect_false(anyNA(res$pop))
  expect_true(all(t(res$pop) >= lo & t(res$pop) <= up))
})
