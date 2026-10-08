## The objective settings must reach every evaluation of the objective: the fit (runDEoptim(), which
## also re-scores the final population) AND the threshold calibration (runSpreadWithoutDEoptim()).
## The calibration used to leave `weighted` at runSpreadWithoutDEoptim()'s default TRUE and never
## passed `sizeLik`, so the threshold came from a different objective than the one it is used in.

recordingRSWD <- function(rec, value) {
  f <- function() {
    rec$args <- as.list(environment())
    value
  }
  formals(f) <- formals(runSpreadWithoutDEoptim)
  f
}

settings <- list(weighted = "sqrt", sizeLik = "t", sizeLikDf = 7, adWeight = 2,
                 yearAreaWeight = "auto", areaDistWeight = "auto", jumpTries = 20, jumpMeanDist = 4)

test_that("the calibration and the fit get the same objective settings", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(c(list(stopIfNoPreRunFit = FALSE, SNLL_FS_thresh = NULL), settings))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, 777))
  sim <- suppressMessages(SpaDES.core::spades(sim))
  for (nm in names(settings)) {
    expect_identical(rec$args[[nm]], settings[[nm]], label = paste("calibration", nm))
    expect_identical(fitRec$deArgs[[nm]], settings[[nm]], label = paste("fit", nm))
  }
  expect_null(rec$args$link)       # logistic3p is the objective's default link, NULL
  expect_null(fitRec$deArgs$link)
})

test_that("runSpreadWithoutDEoptim() hands the settings to the objective", {
  seen <- list()
  local_mocked_bindings(.objfunSpreadFit = function(par, thresh, ...) {
    seen[[length(seen) + 1L]] <<- list(...)
    1000 + thresh
  })
  land <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400)
  fires <- list(year2001 = data.frame(size = c(5L, 12L), date = "year2001", ids = 11:12, cells = c(3L, 56L)),
                year2002 = data.frame(size = c(3L, 10L), date = "year2002", ids = 21:22, cells = c(89L, 12L)))
  ## weighted = "sqrt" used to reach `log(size) ^ weighted` in the rough threshold estimate and error
  suppressMessages(utils::capture.output(do.call(runSpreadWithoutDEoptim, c(list(
    iterThresh = 3L, lower = c(a = 0), upper = c(a = 1), fireSense_spreadFormula = "~ 0 + a",
    flammableRTM = land, annualDTx1000 = list(), nonAnnualDTx1000 = list(), fireBufferedListDT = list(),
    historicalFires = fires, covMinMax = NULL, objfunFireReps = 2L, maxFireSpread = 0.28,
    tests = "SNLL_FS", formulaToFit = "~ 0 + a", mode = "debug", seed = 1L, escapeSizeHa = 50,
    link = "logistic3pUpper"), settings))))
  expect_length(seen, 3L)
  for (nm in names(settings))
    expect_true(all(vapply(seen, function(a) identical(a[[nm]], settings[[nm]]), logical(1))), label = nm)
  expect_true(all(vapply(seen, function(a) identical(a$link, "logistic3pUpper"), logical(1))))
})

test_that("the new terms and jumping are on by default", {
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE, SNLL_FS_thresh = NULL))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, 777))
  sim <- suppressMessages(SpaDES.core::spades(sim))
  for (a in list(rec$args, fitRec$deArgs)) {
    expect_identical(a$yearAreaWeight, "auto")
    expect_identical(a$areaDistWeight, "auto")
    expect_identical(a$jumpTries, 20)
    expect_identical(a$jumpMeanDist, 3)
    expect_identical(a$penaliseRunaways, TRUE)
    expect_identical(a$weighted, FALSE)   # the module's default, now also in the calibration
  }
  for (a in list(rec$args, fitRec$deArgs)) {
    expect_identical(a$runawayEdgeFrac, 0.01)
    expect_identical(a$runawayEdgeMin, 3L)
  }
})

test_that("penaliseRunaways = FALSE reaches the fit, the calibration and the objective", {
  ## FALSE only: the default TRUE is asserted in "the new terms and jumping are on by default", and a
  ## second TRUE call here would hit that test's cached threshold and never reach the mock
  for (value in list(FALSE)) {
    rec <- new.env(); fitRec <- new.env()
    sim <- toySim(list(stopIfNoPreRunFit = FALSE, SNLL_FS_thresh = NULL, penaliseRunaways = value))
    mockFitAndLedger(sim, fitRec)
    mockInModule(sim, runSpreadWithoutDEoptim = recordingRSWD(rec, 777))
    sim <- suppressMessages(SpaDES.core::spades(sim))
    expect_identical(rec$args$penaliseRunaways, value)
    expect_identical(fitRec$deArgs$penaliseRunaways, value)
  }
  expect_identical(formals(runSpreadWithoutDEoptim)$penaliseRunaways, TRUE)
  seen <- list()
  local_mocked_bindings(.objfunSpreadFit = function(par, thresh, ...) {
    seen[[length(seen) + 1L]] <<- list(...)
    1000 + thresh
  })
  land <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400)
  fires <- list(year2001 = data.frame(size = c(5L, 12L), date = "year2001", ids = 11:12, cells = c(3L, 56L)))
  suppressMessages(utils::capture.output(runSpreadWithoutDEoptim(
    iterThresh = 3L, lower = c(a = 0), upper = c(a = 1), fireSense_spreadFormula = "~ 0 + a",
    flammableRTM = land, annualDTx1000 = list(), nonAnnualDTx1000 = list(), fireBufferedListDT = list(),
    historicalFires = fires, covMinMax = NULL, objfunFireReps = 2L, maxFireSpread = 0.28,
    tests = "SNLL_FS", formulaToFit = "~ 0 + a", mode = "debug", seed = 1L, escapeSizeHa = 50,
    penaliseRunaways = FALSE, runawayEdgeFrac = 0.5, runawayEdgeMin = 7L)))
  expect_true(length(seen) > 0)
  expect_true(all(vapply(seen, function(a) identical(a$runawayEdgeFrac, 0.5) && identical(a$runawayEdgeMin, 7L),
                         logical(1))))
  expect_true(all(vapply(seen, function(a) identical(a$penaliseRunaways, FALSE), logical(1))))
})
