## The module's defaults are the decided configuration (trial lp05 of FireSenseTesting's expt.R), so a fit launched
## without a trial name does not fit the October objective. Trials set these explicitly; this catches drift.
test_that("the parameter defaults are the decided configuration", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)$parameters
  def <- stats::setNames(md$default, md$paramName)
  expect_identical(def$spreadLink, "powerLawSmooth")
  expect_identical(def$powerLawShape, 0.5)
  expect_identical(def$lowerAsymptote, fireSenseUtils::spreadProbFloor) # 0.13, as lp05
  expect_identical(def$wholeLandscapeSpread, TRUE)
  expect_identical(def$objfunFireReps, 50L)
  expect_identical(def$gateRefusal, "simulate")
  expect_identical(def$annualAreaYearWeights, "equal")
  expect_identical(def$DEoptimTests, "crps")
  expect_identical(def$crpsWeights, list(fire = 1, adTest = 0.34, areaDist = 2.6, totalArea = 164,
                                              top10 = 49.5, annualArea = 520))
  expect_identical(def$fireSpreadSDBounds, NA_real_)
  expect_identical(def$yearSpreadSDBounds, c(0, 10)) # not settled: unchanged
  expect_identical(def$coefSigns, "free")
})
