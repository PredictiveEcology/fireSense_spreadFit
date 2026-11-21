## The objective is a callee of runDEoptim(); Cache() digests the call's arguments, not what the objective
## does. Its bodies go in the fit's key, so a change in how fits are scored is a cache
## miss (2026-10-01: fireSenseUtils' runaway censoring changed, and the old under-burning fits must not be served).
test_that("the fit's cache key includes the objective's bodies", {
  ob <- objectiveBodies()
  ## the two functions DEoptim evaluates, and what a hand list once held ...
  expect_true(all(c(".objfunSpreadFit", "objFunInner", "crpsEnsemble", "crpsFireSizes", "crpsTotalArea",
                    "crpsTopFires", "crpsWeightsDefault", "spreadProbFromIntegerCovs", "spreadProbCovsOfYear",
                    "spreadProbFromCovs", "logisticAll", "powerLawLink", "powerLawSmoothLink", "richardsLink",
                    "linearLinkPars", "linearLinkProb", "spreadLinkLower", "dataYearBefore",
                    "addBufferEdge", "bufferEdge", "bufferOpen", "sizeLikT",
                    "adStatistic", "adWeightAuto", "pooledFireSizes", "yearAreaNLL", "areaWeightedCvM",
                    "pitMidRank", "adUniform") %in% names(ob)))
  ## ... and what it had missed: the annual-area and tail-share terms, the climate ceiling and activation, the split
  ## climate terms, the link's parameters, the gates, the crop to the cells a fire reaches
  expect_true(all(c("crpsAnnualArea", "crpsTopShare", "spreadLinkInputs", "spreadClimateSplit", "splitSpreadPar",
                    "linkParsMatrix", "unburnableLinkPars", "spreadProbGateTest", "cropToCells") %in% names(ob)))
  expect_identical(names(ob), sort(names(ob))) # one order, so the key is stable
  expect_identical(ob, objectiveBodies())
  ## a function the objective does not reach is not in the key
  expect_false("calibratePowerLaw" %in% names(ob))
  ## which replicates are cut at the buffers, the per-fire score, and the terms added to them
  for (f in c("sizeLikT", "adStatistic", "yearAreaNLL", "areaWeightedCvM", "pooledFireSizes",
              "pitMidRank", "adUniform", "crpsEnsemble", "crpsFireSizes", "crpsTotalArea", "crpsTopFires"))
    expect_identical(ob[[f]], deparse(utils::getFromNamespace(f, "fireSenseUtils")), label = f)
  expect_identical(ob$bufferOpen, deparse(fireSenseUtils::bufferOpen))
  ## the power-law link, and the dispatcher that reaches it
  expect_identical(ob$powerLawLink, deparse(fireSenseUtils::powerLawLink))
  expect_identical(ob$powerLawSmoothLink, deparse(fireSenseUtils::powerLawSmoothLink))
  expect_identical(ob$richardsLink, deparse(fireSenseUtils::richardsLink))
  expect_identical(ob$linearLinkProb, deparse(fireSenseUtils::linearLinkProb))
  expect_identical(ob$logisticAll, deparse(fireSenseUtils::logisticAll))
  ## the spread link itself: a change in how spreadProb is computed from the covariates is a cache miss
  expect_identical(ob$spreadProbFromCovs, deparse(fireSenseUtils::spreadProbFromCovs))
  ## the covariates the objective reads are part of it: a change in which vegetation table a fire year uses
  ## must not serve a fit cached under the old rule
  expect_identical(ob$spreadProbFromIntegerCovs, deparse(fireSenseUtils::spreadProbFromIntegerCovs))
  expect_identical(ob$dataYearBefore, deparse(fireSenseUtils::dataYearBefore))
  expect_identical(ob$objFunInner, deparse(utils::getFromNamespace("objFunInner", "fireSenseUtils")))
  fit <- paste(readLines("../../R/fitSpread.R", warn = FALSE), collapse = "\n")
  expect_match(fit, ".cacheExtra = list(fnName, objectiveBodies())", fixed = TRUE)
})

test_that("changing a function the objective reaches only indirectly changes the key; an unrelated one does not", {
  before <- objectiveBodies()
  edited <- function(f) {
    g <- utils::getFromNamespace(f, "fireSenseUtils")
    body(g) <- call("{", quote(.editedForTest <- TRUE), body(g))
    g
  }
  ## crpsAnnualArea is called by the objective's per-year terms, linkParsMatrix by the link; neither was in the key
  for (f in c("crpsAnnualArea", "linkParsMatrix", "spreadLinkInputs")) {
    after <- local({
      do.call(testthat::local_mocked_bindings, c(stats::setNames(list(edited(f)), f), list(.package = "fireSenseUtils")))
      objectiveBodies()
    })
    expect_false(identical(before, after), label = f)
    expect_identical(names(after), names(before))
  }
  unrelated <- testthat::with_mocked_bindings(objectiveBodies(), calibratePowerLaw = edited("calibratePowerLaw"),
                                              .package = "fireSenseUtils")
  expect_identical(unrelated, before)
})

test_that("the emulator's cache key holds every module function the start population reaches, found, not listed", {
  eb <- emulatorBodies()
  expect_true(all(c("emulatorStartPopulation", "emulatorDefaults", "emulatorControlOf", "emulatorSpace", "emulatorTrustBox",
                    "emulatorAtBoxEdge", "emulatorTrustStep", "emulatorNoiseSD", "emulatorPointKeys", "emulateObjective",
                    "emulatorScale", "fitEmulatorGP", "emulatorConverged", "selectEmulatorBatch",
                    "emulatorSpreadPopulation", "emulatorEvaluator", "objectiveFnArgs") %in% names(eb)))
  expect_identical(names(eb), sort(names(eb)))
  ## a function the emulator does not reach is not in the key
  expect_false("estimateSpreadParams" %in% names(eb))
  edited <- function(f) {
    g <- get(f, mode = "function")
    body(g) <- call("{", quote(.editedForTest <- TRUE), body(g))
    g
  }
  for (f in c("emulatorNoiseSD", "objectiveFnArgs", "emulatorTrustStep")) {
    after <- local({
      do.call(testthat::local_mocked_bindings, stats::setNames(list(edited(f)), f))
      emulatorBodies()
    })
    expect_false(identical(eb, after), label = f)
  }
  unrelated <- testthat::with_mocked_bindings(emulatorBodies(), estimateSpreadParams = edited("estimateSpreadParams"))
  expect_identical(unrelated, eb)
})
