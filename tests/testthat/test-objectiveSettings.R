## The objective settings must reach every evaluation of the objective: the fit (runDEoptim(), which
## also re-scores the final population) AND the debug event (runSpreadWithoutDEoptim()).
## Helpers: debugArgs() and fitArgs() in helper-debugEvent.R.

settings <- list(weighted = "sqrt", sizeLik = "t", sizeLikDf = 7, adWeight = 2,
                 yearAreaWeight = "auto", areaDistWeight = "auto", jumpTries = 20, jumpMeanDist = 4, annualAreaYearWeights = "climate")

test_that("the debug call and the fit get the same objective settings", {
  dbg <- debugArgs(settings)
  fit <- fitArgs(settings)
  for (nm in names(settings)) {
    expect_identical(dbg[[nm]], settings[[nm]], label = paste("debug", nm))
    expect_identical(fit[[nm]], settings[[nm]], label = paste("fit", nm))
  }
  expect_null(dbg$link)       # logistic3p is the objective's default link, NULL
  expect_null(fit$link)
})

test_that("runSpreadWithoutDEoptim() hands the settings to the objective", {
  seen <- list()
  local_mocked_bindings(.objfunSpreadFit = function(par, ...) {
    seen[[length(seen) + 1L]] <<- list(...)
    1000
  })
  land <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400)
  fires <- list(year2001 = data.frame(size = c(5L, 12L), date = "year2001", ids = 11:12, cells = c(3L, 56L)),
                year2002 = data.frame(size = c(3L, 10L), date = "year2002", ids = 21:22, cells = c(89L, 12L)))
  suppressMessages(utils::capture.output(do.call(runSpreadWithoutDEoptim, c(list(
    nSets = 3L, lower = c(a = 0), upper = c(a = 1), fireSense_spreadFormula = "~ 0 + a",
    flammableRTM = land, annualDTx1000 = list(), nonAnnualDTx1000 = list(), fireBufferedListDT = list(),
    historicalFires = fires, covMinMax = NULL, objfunFireReps = 2L, maxFireSpread = 0.28,
    tests = "crps", formulaToFit = "~ 0 + a", seed = 1L, escapeSizeHa = 50, doObjFunAssertions = FALSE,
    link = "logistic3pUpper"), settings))))
  expect_length(seen, 3L)
  for (nm in names(settings))
    expect_true(all(vapply(seen, function(a) identical(a[[nm]], settings[[nm]]), logical(1))), label = nm)
  expect_true(all(vapply(seen, function(a) identical(a$link, "logistic3pUpper"), logical(1))))
})

test_that("the new terms and jumping are on by default", {
  for (a in list(debugArgs(), fitArgs())) {
    expect_identical(a$yearAreaWeight, "auto")
    expect_identical(a$areaDistWeight, "auto")
    expect_identical(a$jumpTries, 20)
    expect_identical(a$jumpMeanDist, 3)
    expect_identical(a$penaliseRunaways, TRUE)
    expect_identical(a$weighted, FALSE)
    expect_identical(a$runawayEdgeFrac, 0.01)
    expect_identical(a$runawayEdgeMin, 3L)
    expect_null(a$runawayBufferMultiple)    # NA: fireSenseUtils does not warn about it
  }
})

test_that("runawayBufferMultiple = NA reaches the fit and the debug call as NULL", {
  dbg <- debugArgs(list(runawayBufferMultiple = NA_real_))
  fit <- fitArgs(list(runawayBufferMultiple = NA_real_))
  expect_true("runawayBufferMultiple" %in% names(dbg))
  expect_null(dbg$runawayBufferMultiple)
  expect_true("runawayBufferMultiple" %in% names(fit))
  expect_null(fit$runawayBufferMultiple)
  expect_null(formals(runSpreadWithoutDEoptim)$runawayBufferMultiple)
})

test_that("penaliseRunaways = FALSE reaches the fit, the debug call and the objective", {
  expect_identical(debugArgs(list(penaliseRunaways = FALSE))$penaliseRunaways, FALSE)
  expect_identical(fitArgs(list(penaliseRunaways = FALSE))$penaliseRunaways, FALSE)
  expect_identical(formals(runSpreadWithoutDEoptim)$penaliseRunaways, TRUE)
  seen <- list()
  local_mocked_bindings(.objfunSpreadFit = function(par, ...) {
    seen[[length(seen) + 1L]] <<- list(...)
    1000
  })
  land <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400)
  fires <- list(year2001 = data.frame(size = c(5L, 12L), date = "year2001", ids = 11:12, cells = c(3L, 56L)))
  suppressMessages(utils::capture.output(runSpreadWithoutDEoptim(
    nSets = 3L, lower = c(a = 0), upper = c(a = 1), fireSense_spreadFormula = "~ 0 + a",
    flammableRTM = land, annualDTx1000 = list(), nonAnnualDTx1000 = list(), fireBufferedListDT = list(),
    historicalFires = fires, covMinMax = NULL, objfunFireReps = 2L, maxFireSpread = 0.28,
    tests = "crps", formulaToFit = "~ 0 + a", seed = 1L, escapeSizeHa = 50, doObjFunAssertions = FALSE,
    penaliseRunaways = FALSE, runawayEdgeFrac = 0.5, runawayEdgeMin = 7L, runawayBufferMultiple = 3)))
  expect_true(length(seen) > 0)
  expect_true(all(vapply(seen, function(a) identical(a$runawayBufferMultiple, 3), logical(1))))
  expect_true(all(vapply(seen, function(a) identical(a$runawayEdgeFrac, 0.5) && identical(a$runawayEdgeMin, 7L),
                         logical(1))))
  expect_true(all(vapply(seen, function(a) identical(a$penaliseRunaways, FALSE), logical(1))))
})

test_that("the CRPS objective's options reach the fit and the debug call, at their defaults too (the decided configuration)", {
  crps <- list(crpsWeights = list(top10 = 5), gateRefusal = "simulate")
  for (a in list(debugArgs(crps), fitArgs(crps))) {
    expect_identical(a$tests, "crps")
    expect_identical(a$crpsWeights, list(top10 = 5))
    expect_identical(a$gateRefusal, "simulate")
  }
  ## the annual-area year weights reach the fit only when not the default, so a fit's cache key is unchanged
  expect_null(fitArgs()$annualAreaYearWeights)
  ## defaults: the decided configuration (test-decidedDefaults.R)
  for (a in list(debugArgs(), fitArgs()))
    expect_identical(a$crpsWeights$annualArea, 520)
})

test_that("the year climate z of the climate year weights and climateSlope reaches the fit and the debug call, standardised, from the ignition covariate CMDsm", {
  ## toy years 2001 and 2002: two years, so z is +/- 1/sqrt(2)
  for (p in list(list(annualAreaYearWeights = "climate"), list(crpsWeights = list(climateSlope = 5)))) {
    for (a in list(debugArgs(p), fitArgs(p))) {
      expect_identical(a$annualClimateCovariate, "CMDsm")
      expect_equal(unname(a$yearClimateZ), c(1, -1) / sqrt(2)) # 2001 mean 23.3 > 2002 mean 22
      expect_named(a$yearClimateZ, c("year2001", "year2002"))
    }
  }
  ## off by default: nothing reaches the fit, so its cache key is unchanged
  expect_null(fitArgs()$yearClimateZ)
  expect_null(fitArgs()$annualClimateCovariate)
})

test_that("a covariate not in the spread data comes from the ignition covariates; with neither, the fit stops with a message", {
  ann <- list(year2001 = data.table::data.table(pixelID = 1:2, CMDsm = c(1, 3)), year2002 = data.table::data.table(pixelID = 1:2, CMDsm = c(5, 7)))
  Par <- list(annualAreaYearWeights = "climate", annualClimateCovariate = "CMDsp", crpsWeights = list(annualArea = 1))
  ign <- data.frame(year = c(2001, 2001, 2002, 2002), CMDsp = c(1, 3, 10, 30), coreFraction = 1)
  r <- yearClimateArgs(Par, ann, ign)
  expect_identical(r$annualClimateCovariate, "CMDsp")
  expect_equal(unname(r$yearClimateZ), c(-1, 1) / sqrt(2))
  expect_error(yearClimateArgs(Par, ann, NULL), "annualClimateCovariate")
  expect_identical(yearClimateArgs(list(annualAreaYearWeights = "equal", crpsWeights = list(climateSlope = 0)), ann, NULL), list())
})
