## The objective is a callee of runDEoptim(); Cache() digests the call's arguments, not what the objective
## does. Its bodies go in the fit's and the calibration's keys, so a change in how fits are scored is a cache
## miss (2026-10-01: fireSenseUtils' runaway censoring changed, and the old under-burning fits must not be served).
test_that("the fit's and the threshold calibration's cache keys include the objective's bodies", {
  ob <- objectiveBodies()
  expect_named(ob, c("objfunSpreadFit", "objFunInner", "crpsEnsemble", "crpsFireSizes", "crpsTotalArea",
                     "crpsTopFires", "crpsWeightsDefault", "spreadProbFromIntegerCovs", "spreadProbCovsOfYear",
                     "spreadProbFromCovs", "logisticAll", "powerLawLink", "dataYearBefore",
                     "addBufferEdge", "bufferEdge", "bufferOpen", "sizeLikT", "kernelLogScore",
                     "adStatistic", "adWeightAuto", "pooledFireSizes", "yearAreaNLL", "areaWeightedCvM",
                     "pitMidRank", "adUniform"))
  ## which replicates are cut at the buffers, the per-fire likelihood and score, and the terms added to them
  for (f in c("kernelLogScore", "sizeLikT", "adStatistic", "yearAreaNLL", "areaWeightedCvM", "pooledFireSizes",
              "pitMidRank", "adUniform", "crpsEnsemble", "crpsFireSizes", "crpsTotalArea", "crpsTopFires"))
    expect_identical(ob[[f]], deparse(utils::getFromNamespace(f, "fireSenseUtils")), label = f)
  expect_identical(ob$bufferOpen, deparse(fireSenseUtils::bufferOpen))
  ## the power-law link, and the dispatcher that reaches it
  expect_identical(ob$powerLawLink, deparse(fireSenseUtils::powerLawLink))
  expect_identical(ob$logisticAll, deparse(fireSenseUtils::logisticAll))
  ## the spread link itself: a change in how spreadProb is computed from the covariates is a cache miss
  expect_identical(ob$spreadProbFromCovs, deparse(fireSenseUtils::spreadProbFromCovs))
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
