## Parameter classes and defaults. test-metadata.R pins the names of parameters, inputs and outputs and
## the classes of inputs and outputs; this pins what a user's script relies on when it does NOT set a
## parameter. A changed default changes every run that leaves it alone.

md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
p <- md$parameters
cls <- stats::setNames(vapply(p$paramClass, paste, character(1), collapse = "|"), p$paramName)
def <- stats::setNames(p$default, p$paramName)

test_that("parameter classes", {
  expect_identical(
    cls[sort(names(cls), method = "radix")],
    (function(x) x[sort(names(x), method = "radix")])(c(.c = "numeric", .plotInterval = "numeric", .plots = "character|logical", .plotSize = "list", .runInitialTime = "numeric",
      .studyAreaName = "character", .useCache = "logical|character",
      cores = "integer", DEoptimControl = "list", DEoptimTests = "character",
      doObjFunAssertions = "logical", heldOutFold = "integer", initialpop = "numeric", iterDEoptim = "integer",
      iterThresh = "integer", thresholdMargin = "numeric", libPathDEoptim = "character", lower = "numeric",
      maxFireSpread = "numeric", mode = "character", mutuallyExclusiveCols = "list",
      nCoresNeeded = "integer", objFunCoresInternal = "integer",
      objfunFireReps = "integer", refitExisting = "logical",
      .rep = "integer", rescaleAll = "logical", SNLL_FS_thresh = "integer",
      spreadFitFilename = "character", spreadFitGoogleDriveFolder = "character",
      stopIfNoPreRunFit = "logical", strategy = "integer", trace = "numeric", upper = "numeric",
      upperAndLowerVal = "numeric", upperAndLowerValFuel = "numeric", useCache_DE = "logical",
      covFixedRange = "list",
      link = "character", profileReps = "integer", simulateMembers = "integer",
      sizeLik = "character", sizeLikDf = "numeric", escapeSizeHa = "numeric", weighted = "logical|character",
      adWeight = "character|numeric", yearAreaWeight = "numeric|character",
      areaDistWeight = "numeric|character", penaliseRunaways = "logical", runawayEdgeFrac = "numeric", runawayEdgeMin = "integer", jumpTries = "numeric", jumpMeanDist = "numeric",
      upperTailBounds = "numeric", yearSpreadSDBounds = "numeric",
      verbose = "numeric", visualizeDEoptim = "Path"))
  )
})

test_that("the defunct DEoptim-retrieval parameters and output are gone", {
  for (nm in c("cacheId_DE", "cloudFolderID_DE", "useCloud_DE", "onlyLoadDEOptim", "NP",
               ".runInterval", ".saveInitialTime", ".saveInterval", "urlDEOptimObject"))
    expect_false(nm %in% p$paramName, label = nm)
  expect_false("fireSense_SpreadFitted" %in% md$outputObjects$objectName)
})

test_that("DEoptim defaults are the tuning study's best and a production fit's settings", {
  ## Tuning study (ELF 13.1, NP 110, 150 generations, 3 seeds each, 2026-09-15): strategy 6 with p = 0.1 had the
  ## lowest median best value and stayed lowest after re-scoring the best members 10 times. Provisional.
  expect_identical(def$strategy, 6L)
  ## c = 0: `iterStep` is hard-coded to 1, so the adaptation restarts every generation, and c > 0 would
  ## describe a run that never happens
  expect_identical(def$.c, 0)
  ## a ceiling: clusters (>= 0.0.46) stops a fit once its population median stops improving
  expect_identical(def$iterDEoptim, 5000L)
  expect_identical(def$objfunFireReps, 50L)
})

test_that("defaults that decide whether and how a fit runs", {
  expect_identical(def$stopIfNoPreRunFit, TRUE)
  expect_identical(def$refitExisting, FALSE)
  expect_identical(def$mode, "fit")
  expect_identical(def$DEoptimTests, c("adTest", "SNLL_FS"))
  expect_identical(def$rescaleAll, TRUE)
  expect_identical(def$doObjFunAssertions, TRUE)
  expect_identical(def$useCache_DE, TRUE)
  expect_identical(def$.useCache, "init")
  expect_identical(def$mutuallyExclusiveCols, list(youngAge = c("class", "nonForest")))
})

test_that("the runaway rule's parameters exist with their defaults", {
  expect_identical(def$runawayEdgeFrac, 0.01)
  expect_identical(def$runawayEdgeMin, 3L)
})

test_that("numeric defaults", {
  expect_identical(def$iterThresh, 96L)
  expect_identical(def$objFunCoresInternal, 1L)
  expect_identical(def$cores, 1L)
  expect_identical(def$.rep, 1L)
  expect_identical(def$trace, 1L)
  expect_identical(def$verbose, 1)
  expect_identical(def$maxFireSpread, 0.28)
  expect_identical(def$upperAndLowerVal, 50)
  expect_identical(def$upperAndLowerValFuel, 100)
  expect_identical(def$covFixedRange, list(CMDsm = c(0, 100), CMD = c(0, 100), CMDsp = c(0, 100), cumMDC = c(0, 100)))
  expect_identical(def$.plotSize, list(height = 1600, width = 2000))
})

test_that("defaults of the objective and of the post-fit diagnostics", {
  ## the likelihood the 2026-09-21 cross-validation chose; the fits before this used "kde" with a
  ## log(size) weight only because the module could not ask for anything else
  expect_identical(def$sizeLik, "t")
  expect_identical(def$sizeLikDf, 5)
  expect_identical(def$weighted, FALSE)
  expect_identical(def$adWeight, "auto")
  ## the annual-area and area-distribution terms and jumping are on by default
  expect_identical(def$yearAreaWeight, "auto")
  expect_identical(def$areaDistWeight, "auto")
  expect_identical(def$jumpTries, 20)
  expect_identical(def$jumpMeanDist, 3)
  expect_identical(def$link, "logistic3p")
  expect_identical(def$upperTailBounds, c(-1, 1))
  expect_identical(def$yearSpreadSDBounds, c(0, 1))   # the per-year random effect is on by default
  expect_identical(def$profileReps, 10L)
  expect_identical(def$.plotInterval, 25)       # DEoptim generations between progress figures
  expect_identical(def$simulateMembers, 10L)
})

test_that("defaults that are 'not set'", {
  for (nm in c(".plots", "initialpop", "nCoresNeeded", "SNLL_FS_thresh"))
    expect_null(def[[nm]], label = nm)
  for (nm in c("lower", "upper"))
    expect_identical(def[[nm]], NA, label = nm)
  ## heldOutFold NA is the "unchanged behaviour" default; see test-initDecision.R
  expect_identical(def$heldOutFold, NA)
})

test_that("the ledger's location", {
  expect_identical(def$spreadFitFilename, "latest")
  expect_identical(def$spreadFitGoogleDriveFolder,
                   "https://drive.google.com/drive/folders/1X9-mRjyLMNpgkP_cfqhbr_AQEPOsVCHf")
})

test_that("the ledger's list-columns, which prediction reads by name", {
  expect_identical(fireSenseUtils::spreadFitAdditionalColNamesTxt,
                   c("numIterations", "objFunVal", "params", "sppEquiv", "nonForestedLCCGroups",
                     "missingLCCgroup", "covMinMax_spread"))
})

test_that("only studyArea has a sourceURL", {
  io <- md$inputObjects
  expect_identical(io$objectName[!is.na(io$sourceURL)], "studyArea")
})
