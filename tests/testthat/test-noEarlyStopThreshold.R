## When every calibration trial fails, the fit must run without the early stop (thresh = Inf), not
## with NA: NA reaches the objective's `SNLL > thresh * 2` and kills every cluster node
## (ELF 6.2.1, heldOutFold 1).

test_that("every trial failed -> Inf, with a message", {
  thresh <- suppressWarnings(pickThreshold(c(1e6, 1e6)))
  expect_message(out <- noEarlyStopThreshold(thresh, runName = "6.2.1_cvFold1", nTrials = 2L),
                 "no SNLL threshold calibrated for 6.2.1_cvFold1: every one of 2 trials failed")
  expect_identical(out, Inf)
})

test_that("a cached NA is not reused as NA", {
  expect_identical(suppressMessages(noEarlyStopThreshold(NA_real_, "x", 96L)), Inf)
  expect_identical(suppressMessages(noEarlyStopThreshold(NA, "x", 96L)), Inf)
})

test_that("a calibrated threshold and NULL (debug mode) pass through silently", {
  expect_silent(expect_identical(noEarlyStopThreshold(600, "x", 96L), 600))
  expect_silent(expect_null(noEarlyStopThreshold(NULL, "x", 96L)))
})

test_that("the consumer applies it after the Cache, and fitSpread() refuses an NA threshold", {
  src <- paste(readLines("../../fireSense_spreadFit.R", warn = FALSE), collapse = "\n")
  expect_match(src, "noEarlyStopThreshold\\(thresh, runName")
  expect_error(fitSpread(sim = NULL, covs = NULL, thresh = NA_real_, runName = "x"), "`thresh` is NA")
})
