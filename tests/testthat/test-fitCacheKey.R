## The objective is a callee of runDEoptim(); Cache() digests the call's arguments, not what the objective
## does. Its bodies go in the fit's and the calibration's keys, so a change in how fits are scored is a cache
## miss (2026-10-01: fireSenseUtils' runaway censoring changed, and the old under-burning fits must not be served).
test_that("the fit's and the threshold calibration's cache keys include the objective's bodies", {
  ob <- objectiveBodies()
  expect_named(ob, c("objfunSpreadFit", "objFunInner", "spreadProbFromIntegerCovs", "spreadProbCovsOfYear",
                     "dataYearBefore"))
  ## the covariates the objective reads are part of it: a change in which vegetation table a fire year uses
  ## must not serve a fit or threshold (possibly NA) cached under the old rule
  expect_identical(ob$spreadProbFromIntegerCovs, deparse(fireSenseUtils::spreadProbFromIntegerCovs))
  expect_identical(ob$dataYearBefore, deparse(fireSenseUtils::dataYearBefore))
  expect_identical(ob$objFunInner, deparse(utils::getFromNamespace("objFunInner", "fireSenseUtils")))
  fit <- paste(readLines("../../R/fitSpread.R", warn = FALSE), collapse = "\n")
  expect_match(fit, ".cacheExtra = list(fnName, objectiveBodies())", fixed = TRUE)
  thr <- paste(readLines("../../fireSense_spreadFit.R", warn = FALSE), collapse = "\n")
  expect_match(thr, "objective = objectiveBodies()", fixed = TRUE)
})
