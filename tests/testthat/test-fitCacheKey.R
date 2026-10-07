## The objective is a callee of runDEoptim(); Cache() digests the call's arguments, not what the objective
## does. Its bodies go in the fit's and the calibration's keys, so a change in how fits are scored is a cache
## miss (2026-10-01: fireSenseUtils' runaway censoring changed, and the old under-burning fits must not be served).
test_that("the fit's and the threshold calibration's cache keys include the objective's bodies", {
  ob <- objectiveBodies()
  expect_named(ob, c("objfunSpreadFit", "objFunInner"))
  expect_identical(ob$objFunInner, deparse(utils::getFromNamespace("objFunInner", "fireSenseUtils")))
  fit <- paste(readLines("../../R/fitSpread.R", warn = FALSE), collapse = "\n")
  expect_match(fit, ".cacheExtra = list(fnName, objectiveBodies())", fixed = TRUE)
  thr <- paste(readLines("../../fireSense_spreadFit.R", warn = FALSE), collapse = "\n")
  expect_match(thr, "objective = objectiveBodies()", fixed = TRUE)
})
