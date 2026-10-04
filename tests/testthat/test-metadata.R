## The module's metadata is its public contract: a project using this module binds
## to these object names and classes. Renaming or retyping one breaks every caller,
## which is exactly the class of change the raster -> terra migration makes, so it is
## worth asserting here rather than discovering downstream.
##
## When a change is deliberate, update this file in the same commit and bump the
## module version to match: removed, renamed or retyped is a MAJOR bump.

test_that("module metadata parses", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_type(md, "list")
  expect_identical(md$name, moduleName)
})

test_that("inputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  inputs <- stats::setNames(md$inputObjects$objectClass, md$inputObjects$objectName)
  expect_identical(
    inputs[order(names(inputs))],
    c(.ELFind                                = "character",
      .runName                               = "character",
      fireBufferedListDT                     = "list",
      fireSense_annualSpreadFitCovariates    = "data.table",
      fireSense_nonAnnualSpreadFitCovariates = "data.table",
      fireSense_spreadFormula                = "character",
      parsKnown                              = "numeric",
      rasterToMatch                          = "SpatRaster",
      spreadFirePoints                       = "sf",
      spreadFitAdditionalColNames            = "character",
      studyArea                              = "sf")
  )
})

test_that("outputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  outputs <- stats::setNames(md$outputObjects$objectClass, md$outputObjects$objectName)
  expect_identical(
    outputs[order(names(outputs))],
    c(covMinMax_spread          = "data.table",
      DE                        = "data.table",
      fsSpreadFit_hists         = "ggplot",
      lociList                  = "list",
      spreadFitConvergence      = "data.table",
      spreadFitHeldOut          = "list",
      spreadFitIdentifiability  = "data.table",
      spreadFitLinkSaturation   = "data.table",
      spreadFitProfile          = "data.table",
      spreadFitRescore          = "data.table",
      spreadFitSizes            = "data.table",
      studyAreaWithSpreadParams = "sf")
  )
})

test_that("parameters are the expected names", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(
    sort(md$parameters$paramName),
    sort(c(".c", ".plotInterval", ".plots", ".plotSize", ".runInitialTime", ".studyAreaName",
           ".useCache", "cores", "covFixedRange", "DEoptimControl", "DEoptimTests", "doObjFunAssertions",
           "heldOutFold", "initialpop", "iterDEoptim", "iterThresh", "thresholdMargin", "libPathDEoptim",
           "link", "lower", "maxFireSpread", "mode", "mutuallyExclusiveCols", "nCoresNeeded",
           "objFunCoresInternal", "objfunFireReps", ".rep", "adWeight", "profileReps",
           "simulateMembers", "sizeLik", "sizeLikDf", "escapeSizeHa", "yearAreaWeight", "areaDistWeight", "penaliseRunaways", "runawayEdgeFrac", "runawayEdgeMin",
           "jumpTries", "jumpMeanDist", "upperTailBounds", "weighted", "yearSpreadSDBounds",
           "refitExisting", "rescaleAll", "SNLL_FS_thresh", "spreadFitFilename",
           "spreadFitGoogleDriveFolder", "stopIfNoPreRunFit", "strategy", "trace",
           "upper", "upperAndLowerVal", "upperAndLowerValFuel", "useCache_DE",
           "verbose", "visualizeDEoptim"))
  )
})

test_that("the required clusters has the fixes a fit on the fleet needs", {
  ## clusters 0.0.41 (2026-09-15): clusterSetup() no longer sends the cluster object to every worker (a fit
  ## stalled for hours), works without reproducible attached or ~/.ssh/config, and picks tunnel ports below
  ## the ephemeral range (a 110-worker build hung). With a lower floor, Require keeps an installed clusters
  ## that has none of these (the fleet had 0.0.31).
  ## clusters 0.0.42: DEoptimIterative2() runs DEoptim with c = 0. With c > 0, DEoptim's F adaptation turns
  ## every trial vector into NaN once a call's first generation has no success.
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  clustersReq <- grep("/clusters@", unlist(md$reqdPkgs), value = TRUE)
  expect_length(clustersReq, 1L)
  ## clusters 0.0.46: DEoptimIterative2() stops a fit when the population median has stopped improving, so
  ## `iterDEoptim` (5000) is a ceiling, not the run length.
  expect_true(package_version(sub(".*>=\\s*([0-9.]+).*", "\\1", clustersReq)) >= "0.0.46")
  ## clusters 0.0.52: DEoptimIterative2(plotEvery), which `.plotInterval` reaches through runDEoptim()
  expect_true(package_version(sub(".*>=\\s*([0-9.]+).*", "\\1", clustersReq)) >= "0.0.52")
})
