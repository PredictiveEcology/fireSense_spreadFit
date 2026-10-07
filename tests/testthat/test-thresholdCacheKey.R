## estimateSNLLThresholdPostLargeFires() cached runSpreadWithoutDEoptim() with
## omitArgs = c("objfunFireReps", "mode"), but both change the result: `mode` picks a
## branch that returns a different type (a numeric threshold vs NULL), and
## `objfunFireReps` changes the threshold's value. Leaving them out of the key served a
## debug-mode result to a fit-mode run, and a 5-replicate threshold to a 25-replicate run.
##
## The call needs a full simList, so this checks the function's source instead of running
## it: parsed, not sourced, so no module machinery is needed.
test_that("the SNLL threshold cache key does not omit mode or objfunFireReps", {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_spreadFit.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name("estimateSNLLThresholdPostLargeFires")), exprs)
  expect_length(def, 1L)
  code <- paste(deparse(def[[1]][[3]]), collapse = "\n")
  expect_match(code, "runSpreadWithoutDEoptim", fixed = TRUE)
  expect_false(grepl("omitArgs\\s*=\\s*c\\([^)]*\"(mode|objfunFireReps)\"", code))
})

## The runaway rule only changes how quickly DEoptim moves away from an unlucky draw, not what a fit
## means, so fits cached under the 1-cell rule must not rerun: both parameters are omitted from the key.
test_that("the DEoptim fit's cache key omits runawayEdgeFrac and runawayEdgeMin", {
  code <- paste(readLines(testthat::test_path("..", "..", "R", "fitSpread.R")), collapse = "\n")
  expect_match(code, "omitArgs\\s*=\\s*c\\([^)]*\"runawayEdgeFrac\"[^)]*\"runawayEdgeMin\"")
})

## The calibration is passed the rule too (the same objective as the fit) but, like the fit, omits it
## from its key: thresholds cached under the 1-cell rule stay valid.
test_that("the threshold calibration passes runawayEdgeFrac and runawayEdgeMin and omits them from its key", {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_spreadFit.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name("estimateSNLLThresholdPostLargeFires")), exprs)
  code <- paste(deparse(def[[1]][[3]]), collapse = "\n")
  expect_match(code, "runawayEdgeFrac = P(sim)$runawayEdgeFrac", fixed = TRUE)
  expect_match(code, "runawayEdgeMin = P(sim)$runawayEdgeMin", fixed = TRUE)
  expect_match(code, "omitArgs\\s*=\\s*c\\(\"runawayEdgeFrac\",\\s*\"runawayEdgeMin\"[,)]")  # deparse may wrap the line
})

## runawayBufferMultiple changes what a runaway scores, so it is in both keys; NA (NULL to the objective)
## is omitted, so fits and thresholds cached before the parameter existed keep their keys.
test_that("the fit and the calibration pass runawayBufferMultiple, omitted from their keys only when NULL", {
  fit <- paste(readLines(testthat::test_path("..", "..", "R", "fitSpread.R")), collapse = "\n")
  exprs <- parse(testthat::test_path("..", "..", "fireSense_spreadFit.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name("estimateSNLLThresholdPostLargeFires")), exprs)
  thr <- paste(deparse(def[[1]][[3]]), collapse = "\n")
  for (code in list(fit, thr)) {
    expect_match(code, "runawayMultiple <- runawayBufferMultipleOrNULL(P(sim)$runawayBufferMultiple)", fixed = TRUE)
    expect_match(code, "runawayBufferMultiple = runawayMultiple,", fixed = TRUE)
    expect_match(code, "omitNullArgs\\([^)]*runawayBufferMultiple = runawayMultiple\\)")
    expect_false(grepl("omitArgs\\s*=\\s*c\\([^)]*\"runawayBufferMultiple\"", code))
  }
  expect_null(runawayBufferMultipleOrNULL(NA_real_))
  expect_identical(runawayBufferMultipleOrNULL(2), 2)
})
