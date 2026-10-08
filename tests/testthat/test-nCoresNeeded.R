## The DEoptim cluster's size is also its population size: clusters:::.clusterNP() sets NP to the
## number of workers the cluster was built with, discarding any NP the caller asked for. So the only
## way to choose NP is to choose the cluster size, which is `nCoresNeeded` in
## fireSenseUtils::runDEoptim() -- and this module never passed it, so the package default
## (10 * length(lower), i.e. 120 for a 12-parameter fit) always won.
##
## FireSense phase 2 (2026-09-16): measured per-ELF, a generation costs the SLOWEST of NP evaluations,
## and the expected maximum barely moves with NP (68.0 s at NP 120 vs 64.9 s at NP 60 for ELF 4.1).
## Throughput therefore comes from running more ELFs at once, not from shorter generations: at 704
## cores, NP 60 allows 11 concurrent fits where NP 120 allows 5. Choosing NP needs it to be a
## parameter of this module.

test_that("the module declares an nCoresNeeded parameter", {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_spreadFit.R"), keep.source = FALSE)
  code <- paste(vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1)),
                collapse = "\n")
  expect_match(code, "defineParameter\\(\"nCoresNeeded\"")
})

test_that("the runDEoptim() call passes nCoresNeeded from the parameter", {
  ## Parsed, not run: the call needs a full simList.
  exprs <- moduleSource()
  calls <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("runDEoptim"))) calls[[length(calls) + 1L]] <<- x
      for (i in seq_along(x)[-1L]) if (is.call(x[[i]])) walk(x[[i]])
    }
  }
  for (e in exprs) walk(e)
  expect_gte(length(calls), 1L)
  for (cl in calls)
    expect_identical(deparse(cl$nCoresNeeded), "P(sim)$nCoresNeeded")
})
