## `escapeSizeHa` (default 50): the spread model is fitted to escaped fires only, and each simulated fire
## starts as an escaped fire. It must reach the fit AND the threshold calibration, or the calibrated
## threshold comes from a different objective than the one it is used in.

test_that("escapeSizeHa is a numeric parameter, 50 by default", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  p <- md$parameters
  i <- which(p$paramName == "escapeSizeHa")
  expect_length(i, 1L)
  expect_identical(p$default[[i]], 50)
})

test_that("the fit and both threshold calibrations pass escapeSizeHa", {
  exprs <- moduleSource()
  calls <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("runDEoptim")) || identical(x[[1]], as.name("runSpreadWithoutDEoptim")))
        calls[[length(calls) + 1L]] <<- x
      lapply(as.list(x), walk)
    }
    invisible(NULL)
  }
  invisible(lapply(exprs, walk))
  ## the fit, the `debug` event and the `estimateThreshold` event (a packaged module may list a call twice)
  fns <- vapply(calls, function(cl) deparse(cl[[1]]), character(1))
  expect_true(all(c("runDEoptim", "runSpreadWithoutDEoptim") %in% fns))
  expect_gte(sum(fns == "runSpreadWithoutDEoptim"), 2L)
  for (cl in calls)
    expect_identical(deparse(as.list(cl)$escapeSizeHa), "escapeSizeHaOrNULL(P(sim)$escapeSizeHa)")
})

test_that("NULL or NA means the old fit", {
  expect_null(escapeSizeHaOrNULL(NULL))
  expect_null(escapeSizeHaOrNULL(NA))
  expect_identical(escapeSizeHaOrNULL(50), 50)
})

test_that("the threshold calibration hands escapeSizeHa to the objective", {
  seen <- list()
  local_mocked_bindings(.objfunSpreadFit = function(par, thresh, ...) {
    seen[[length(seen) + 1L]] <<- list(...)
    1000 + thresh
  })
  land <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400) # 5.76-ha pixels
  fires <- list(year2001 = data.frame(size = c(5L, 12L), date = "year2001", ids = 11:12, cells = c(3L, 56L)),
                year2002 = data.frame(size = c(3L, 10L), date = "year2002", ids = 21:22, cells = c(89L, 12L)))
  suppressMessages(utils::capture.output(runSpreadWithoutDEoptim(
    iterThresh = 3L, lower = c(a = 0), upper = c(a = 1), fireSense_spreadFormula = "~ 0 + a",
    flammableRTM = land, annualDTx1000 = list(), nonAnnualDTx1000 = list(), fireBufferedListDT = list(),
    historicalFires = fires, covMinMax = NULL, objfunFireReps = 2L, maxFireSpread = 0.28,
    tests = "SNLL_FS", formulaToFit = "~ 0 + a", mode = "debug", seed = 1L, escapeSizeHa = 50)))
  expect_length(seen, 3L)
  expect_true(all(vapply(seen, function(a) identical(a$escapeSizeHa, 50), logical(1))))
})
