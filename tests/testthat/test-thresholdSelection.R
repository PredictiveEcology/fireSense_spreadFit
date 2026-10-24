## Choosing the calibrated early-stop threshold from the trials' own first-block values.
##
## FireSense, 2026-09-30: each trial used to be paired with an independent random threshold and the
## smallest passing one was kept, so the result was luck (ELF 6.2.1 fold 1: non-saturated trials scored
## 201, 212, 226 but the paired thresholds were 17, 209, 144 -> NA; 14.3 fold 1: 533/541 against a
## largest possible threshold of 381 -> NA). Now: threshold = margin x the best usable trial's
## first-block value per fire, the unit fireSenseUtils' objective compares `thresh` with: it bails when
## the first block's per-fire term / its number of fires > thresh (R/objFunSpread.R), that ratio being
## the objective's `firstBlockPerFire`. Each fire's score is about 1-3, so the threshold is not rounded.

test_that("the threshold is 2 x the best usable value; saturated and failed trials are ignored", {
  perFire <- c(2.09, 2.4, 1e6, 0.5, NA, 3.1)
  saturated <- c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)    # 0.5 saturated; 1e6 failed
  expect_identical(pickThreshold(perFire, saturated), 4.18)
  expect_identical(pickThreshold(perFire, saturated, margin = 3), 3 * 2.09)
})

test_that("the threshold is not rounded: per fire, the values are about 1-3", {
  expect_identical(pickThreshold(c(1.15, 50), margin = 2), 2.3)
})

test_that("no usable trial -> NA with a warning, which noEarlyStopThreshold() turns into Inf", {
  expect_warning(out <- pickThreshold(c(1e6, 2e6, NA), c(FALSE, FALSE, FALSE)), "no threshold calibrated")
  expect_identical(out, NA_real_)
  expect_warning(pickThreshold(c(300, 400), c(TRUE, TRUE)), "no threshold calibrated")
  expect_message(res <- noEarlyStopThreshold(out, runName = "x", nTrials = 3L), "every one of 3 trials failed")
  expect_identical(res, Inf)
})

test_that("failVal is the cutoff for a usable value, and is exclusive", {
  expect_identical(pickThreshold(c(99, 100), failVal = 100), 198)
})

test_that("the best usable trial does not bail under the derived threshold (the objective's own comparison)", {
  ## the comparison in fireSenseUtils::.objfunSpreadFit(), first block of nFires fires: its summed term per fire
  bails <- function(blockValue, thresh, nFires) blockValue / nFires > thresh
  for (nFires in c(1L, 7L, 300L)) for (a in c(0.3, 1.15, 2.7, 60)) {
    thresh <- pickThreshold(c(a, 3 * a, 1e6), margin = 2)
    expect_false(bails(a * nFires, thresh, nFires))
    expect_true(bails(2.01 * thresh * nFires, thresh, nFires))   # a worse set still bails
  }
})

test_that("the cache key of the calibration includes the rule's code", {
  src <- paste(readLines("../../fireSense_spreadFit.R", warn = FALSE), collapse = "\n")
  expect_match(src, "\\.cacheExtra = list\\(pickThreshold = deparse\\(pickThreshold\\)")
  expect_match(src, "trialFirstBlock = deparse\\(trialFirstBlock\\)")
})

test_that("a trial's first block is read from the objective's returnTerms fields", {
  seen <- NULL
  fake <- function(par, thresh, returnTerms, ...) {
    seen <<- list(thresh = thresh, returnTerms = returnTerms, dots = list(...))
    c(objective = 1533, firstBlockPerFire = 533, bailed = 1)
  }
  testthat::local_mocked_bindings(.objfunSpreadFit = fake)
  expect_identical(trialFirstBlock(c(a = 1), penaliseRunaways = FALSE),
                   list(perFire = 533, saturated = TRUE))
  expect_identical(seen$thresh, Inf)
  expect_true(seen$returnTerms)
  expect_identical(seen$dots, list(penaliseRunaways = FALSE))
  testthat::local_mocked_bindings(.objfunSpreadFit = function(...) c(objective = 1, firstBlockPerFire = 80, bailed = 0))
  expect_identical(trialFirstBlock(c(a = 1)), list(perFire = 80, saturated = FALSE))
})

test_that("a crashed trial, or an objective without the fields, counts as failed", {
  testthat::local_mocked_bindings(.objfunSpreadFit = function(...) stop("boom"))
  expect_identical(trialFirstBlock(c(a = 1)), list(perFire = NA_real_, saturated = NA, error = "boom"))
  testthat::local_mocked_bindings(.objfunSpreadFit = function(...) 1234)
  expect_match(trialFirstBlock(c(a = 1))$error, "firstBlockPerFire")
})

test_that("when every trial errored the calibration stops with the first error, it does not return NA", {
  errored <- list(trialFirstBlock = list(perFire = NA_real_, saturated = NA, error = "no covariates"),
                  other = list(perFire = NA_real_, saturated = NA, error = "second"))
  expect_error(stopIfAllTrialsErrored(errored), "every calibration trial failed.*no covariates")
  ## one usable trial, or trials that merely bailed (no error), are not an error
  expect_null(stopIfAllTrialsErrored(c(list(list(perFire = 80, saturated = FALSE)), errored)))
  expect_null(stopIfAllTrialsErrored(list(list(perFire = 1e6, saturated = TRUE))))
})

test_that("no printed output is parsed anywhere in the sources", {
  for (f in list.files("../../R", pattern = "\\.R$", full.names = TRUE)) {
    code <- grep("^\\s*#", readLines(f, warn = FALSE), value = TRUE, invert = TRUE)
    expect_false(any(grepl("capture.output|Avg annual|Decent in 1st|Fail! Bailing", code)), info = f)
  }
})
