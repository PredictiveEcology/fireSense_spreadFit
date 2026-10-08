## estimateSpreadParams(): the default DEoptim bounds.
##
## hillSlope1 (the logistic slope) and inflectionPoint1 (its Richards exponent) are fixed at 1, not
## fitted (see NEWS and fireSenseUtils::fixLogisticPars()), so neither is among the default bounds.

form <- "~ 0 + CMDsm + youngAge + class1 + nf"
annual <- list(year2001 = data.table::data.table(pixelID = 1L, CMDsm = 1, youngAge = 0))

test_that("upper bounds: +B for covariates, 0 for youngAge, maxAsymptote first, no hillSlope1 or inflectionPoint1", {
  expect_identical(
    estimateSpreadParams(form, annual, whichBound = "upper", upperAndLower = 9),
    c(maxAsymptote = fireSenseUtils::spreadProbCeiling,
      CMDsm = 9, youngAge = 0, class1 = 9, nf = 9))
  expect_false(any(names(fireSenseUtils::fixedLogisticPars) %in%
                     names(estimateSpreadParams(form, annual, "upper", upperAndLower = 9))))
  expect_false("inflectionPoint1" %in% names(estimateSpreadParams(form, annual, "upper", upperAndLower = 9)))
})

test_that("lower bounds: -B for non-drought covariates and youngAge, 0 for drought-index terms, no hillSlope1 or inflectionPoint1", {
  expect_identical(
    estimateSpreadParams(form, annual, whichBound = "lower", upperAndLower = 9),
    c(maxAsymptote = 0.25,
      CMDsm = 0, youngAge = -9, class1 = -9, nf = -9))
  expect_false("inflectionPoint1" %in% names(estimateSpreadParams(form, annual, "lower", upperAndLower = 9)))
})

test_that("the bound scales with upperAndLower and every lower bound is <= its upper bound", {
  up <- estimateSpreadParams(form, annual, "upper", upperAndLower = 2.5)
  lo <- estimateSpreadParams(form, annual, "lower", upperAndLower = 2.5)
  expect_identical(unname(up[c("CMDsm", "class1", "nf")]), c(2.5, 2.5, 2.5))
  expect_identical(unname(lo[c("youngAge", "class1", "nf")]), c(-2.5, -2.5, -2.5))
  expect_identical(names(up), names(lo))
  expect_true(all(lo <= up))
})

test_that("a term is bounded by name, not by annual-table membership", {
  ## a non-drought term that IS in the annual table stays symmetric, not floored at 0
  form2 <- "~ 0 + CMD + cumMDC + PPT_sm + youngAge + fuel1 + nf"
  annual2 <- list(year2001 = data.table::data.table(pixelID = 1L, CMD = 1, cumMDC = 1,
                                                     PPT_sm = 1, youngAge = 0))
  up <- estimateSpreadParams(form2, annual2, "upper", upperAndLower = 50,
                             fuelTerms = "fuel1", upperAndLowerFuel = 100)
  lo <- estimateSpreadParams(form2, annual2, "lower", upperAndLower = 50,
                             fuelTerms = "fuel1", upperAndLowerFuel = 100)
  expect_identical(unname(c(lo["CMD"], up["CMD"])), c(0, 50))
  expect_identical(unname(c(lo["cumMDC"], up["cumMDC"])), c(0, 50))
  expect_identical(unname(c(lo["PPT_sm"], up["PPT_sm"])), c(-50, 50))
  expect_identical(unname(c(lo["youngAge"], up["youngAge"])), c(-50, 0))
  expect_identical(unname(c(lo["fuel1"], up["fuel1"])), c(0, 100))  # fuel biomass: lower 0, upper unchanged
  expect_identical(unname(c(lo["nf"], up["nf"])), c(-50, 50))

  ## youngAge absent from the annual table: still bounded as youngAge, by name
  noYA <- list(year2001 = data.table::data.table(pixelID = 1L, CMDsm = 1))
  expect_identical(unname(estimateSpreadParams(form, noYA, "upper", 9)["youngAge"]), 0)
  expect_identical(unname(estimateSpreadParams(form, noYA, "lower", 9)["youngAge"]), -9)
  ## no drought or youngAge term in the formula at all
  expect_identical(estimateSpreadParams("~ 0 + class1", annual, "lower", 3),
                   c(maxAsymptote = 0.25, class1 = -3))
})

test_that("term order follows the formula", {
  expect_identical(names(estimateSpreadParams("~ 0 + nf + CMDsm", annual, "upper", 9)),
                   c("maxAsymptote", "nf", "CMDsm"))
})

test_that("whichBound must be 'upper' or 'lower'", {
  expect_error(estimateSpreadParams(form, annual, whichBound = "both", upperAndLower = 9),
               "whichBound %in% c(\"upper\", \"lower\") is not TRUE", fixed = TRUE)
})

test_that("the per-year random effect's sd is bounded last, after an upper-tail term", {
  ## fireSenseUtils::runDEoptim() finds yearSpreadSD by name and requires it last; the objective takes
  ## the logistic parameters by position, so upperTail1 is 2nd (the fixed hillSlope1 and inflectionPoint1 follow maxAsymptote)
  up <- estimateSpreadParams(form, annual, "upper", 9, upperTailBounds = c(-1, 1), yearSpreadSDBounds = c(0, 1))
  lo <- estimateSpreadParams(form, annual, "lower", 9, upperTailBounds = c(-1, 1), yearSpreadSDBounds = c(0, 1))
  expect_identical(names(up)[2], "upperTail1")
  expect_identical(names(up)[length(up)], "yearSpreadSD")
  expect_identical(unname(c(lo["yearSpreadSD"], up["yearSpreadSD"])), c(0, 1))
  expect_identical(names(up), names(lo))
  ## off: no yearSpreadSD at all
  expect_false("yearSpreadSD" %in% names(estimateSpreadParams(form, annual, "upper", 9)))
})
