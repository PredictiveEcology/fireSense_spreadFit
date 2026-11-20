## powerLawPcAtP50 = TRUE fixes p_c (maxAsymptote) at the calibrated p50: both its bounds are p50. Helpers: helper-spreadLink.R.

test_that("maxAsymptoteBoundsOf(atP50 = TRUE) gives p50 for both bounds; FALSE is as before", {
  expect_identical(maxAsymptoteBoundsOf("powerLawSmooth", p50 = 0.26, atP50 = TRUE), c(0.26, 0.26))
  expect_identical(maxAsymptoteBoundsOf("powerLaw", p50 = 0.26, atP50 = TRUE), c(0.26, 0.26))
  expect_identical(maxAsymptoteBoundsOf("powerLaw", p50 = 0.26, atP50 = FALSE), maxAsymptoteBoundsOf("powerLaw", p50 = 0.26))
  expect_error(maxAsymptoteBoundsOf("powerLaw", p50 = 0.5, atP50 = TRUE), "outside")
})

test_that("powerLawPcAtP50 is off by default and, when TRUE, gives the fit p50 as p_c's lower and upper bound", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_false(md$parameters$default[[which(md$parameters$paramName == "powerLawPcAtP50")]])
  out <- fittedCO(list(spreadLink = "powerLawSmooth", powerLawShape = 1, powerLawPcAtP50 = TRUE))
  p50 <- out$sim$powerLawCalibration$p50
  expect_identical(unname(out$rec$deArgs$lower[["maxAsymptote"]]), p50)
  expect_identical(unname(out$rec$deArgs$upper[["maxAsymptote"]]), p50)
  expect_identical(unname(out$rec$deArgs$upper[["CMDsm"]]) > 0, TRUE)   # the other bounds are unchanged
  off <- fittedCO(list(spreadLink = "powerLawSmooth", powerLawShape = 1))
  expect_lt(unname(off$rec$deArgs$lower[["maxAsymptote"]]), p50)
})

test_that("powerLawPcAtP50 = TRUE stops with the ceiling and with a link without p_c50", {
  s <- toySim(list(spreadLink = "powerLaw", spreadClimate = "ceiling", powerLawPcAtP50 = TRUE), toyObjectsCO())
  mockFitAndLedger(s)
  expect_error(suppressMessages(SpaDES.core::spades(s)), "powerLawPcAtP50 = TRUE cannot be combined with spreadClimate")
  s <- toySim(list(spreadLink = "richards", powerLawPcAtP50 = TRUE), toyObjectsCO())
  mockFitAndLedger(s)
  expect_error(suppressMessages(SpaDES.core::spades(s)), "powerLawPcAtP50 = TRUE needs spreadLink")
})

test_that("DEoptim accepts equal lower and upper bounds and keeps the parameter there", {
  skip_if_not_installed("DEoptim")
  f <- function(x) sum((x[2:3] - c(1, -1))^2) + x[1]
  set.seed(1)
  de <- DEoptim::DEoptim(f, lower = c(0.26, -3, -3), upper = c(0.26, 3, 3),
                         control = DEoptim::DEoptim.control(NP = 30, itermax = 15, trace = FALSE))
  expect_true(all(de$member$pop[, 1] == 0.26))
  expect_identical(unname(de$optim$bestmem[1]), 0.26)
})
