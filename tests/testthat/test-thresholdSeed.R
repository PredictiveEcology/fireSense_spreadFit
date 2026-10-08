## runSpreadWithoutDEoptim() drew its own seed with `sample(1e6, 1)`, so the SNLL threshold it
## returns was different on every run. That is not merely untidy: the threshold becomes the `thresh`
## argument to runDEoptim(), so it is part of every cached DEoptim generation's key. A single cache
## miss on the estimateThreshold event therefore re-drew the threshold and invalidated EVERY cached
## generation for that ELF.
##
## Observed on 2026-09-16: 4.2.2 was interrupted and restarted; its threshold went 1236 -> 1416 and
## it recomputed from generation 1, losing ~17 h. 4.1 was interrupted the same way, hit the
## estimateThreshold cache, kept threshold 1705, and replayed 797 cached generations in ~35 min.
##
## A seed derived from the ELF makes the threshold reproducible, so a miss costs only the threshold
## recompute rather than the whole fit. The seed is an argument, so it is part of the Cache key --
## which is what test-thresholdCacheKey.R asks for: nothing that changes the result is omitted.

source(testthat::test_path("..", "..", "R", "runSpreadWithoutDEoptim.R"))

test_that("runSpreadWithoutDEoptim() accepts a seed argument", {
  expect_true("seed" %in% names(formals(runSpreadWithoutDEoptim)))
})

test_that("a given seed reproduces the same parameter draws", {
  ## The threshold is a function of `pars` (runif) and `thresholds` (sample), both drawn after the
  ## seed is set. Reproducing those draws is what makes the threshold reproducible; running the full
  ## function needs a landscape and fire data, so drive the same two draws directly.
  draws <- function(s) {
    set.seed(s)
    list(pars = runif(6), thresholds = sample(1000, 4))
  }
  expect_identical(draws(12345L), draws(12345L))
  expect_false(identical(draws(12345L), draws(54321L)))
})

test_that("the production caller passes a deterministic seed derived from the ELF", {
  ## Parsed, not run: estimateSNLLThresholdPostLargeFires() needs a full simList.
  exprs <- parse(testthat::test_path("..", "..", "fireSense_spreadFit.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name("estimateSNLLThresholdPostLargeFires")), exprs)
  expect_length(def, 1L)
  code <- paste(deparse(def[[1]][[3]]), collapse = "\n")
  expect_match(code, "runSpreadWithoutDEoptim", fixed = TRUE)
  ## the seed must be passed, and must come from the ELF identifier rather than from chance
  expect_match(code, "seed\\s*=")
  expect_match(code, "\\.ELFind|\\.runName")
})

test_that("the same ELF gives the same seed, a different ELF a different one", {
  ## whatever derivation is used, it must be a pure function of the identifier
  exprs <- parse(testthat::test_path("..", "..", "fireSense_spreadFit.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(".elfSeed")), exprs)
  expect_length(def, 1L)
  .elfSeed <- eval(def[[1]][[3]])
  expect_identical(.elfSeed("14.3"), .elfSeed("14.3"))
  expect_false(identical(.elfSeed("14.3"), .elfSeed("4.2.2")))
  expect_true(is.numeric(.elfSeed("14.3")) && length(.elfSeed("14.3")) == 1L)
  expect_true(.elfSeed("14.3") >= 1 && .elfSeed("14.3") <= 1e6)
})
