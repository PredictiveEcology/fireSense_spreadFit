## DEoptim settings reach DEoptim. Eliot, 2026-09-15: "Any user passed args should pass into the
## DEoptim processes." The module passed strategy, trace, initialpop and .c only, and .c did not
## reach DEoptim at all; `DEoptimControl` carries any other DEoptim.control() setting (CR, F, p,
## reltol, ...) through fireSenseUtils::runDEoptim() to clusters::clusterSetup().

test_that("DEoptimControl is a list parameter; by default it sets only p, for strategy 6", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  p <- md$parameters
  i <- which(p$paramName == "DEoptimControl")
  expect_length(i, 1L)
  expect_identical(p$paramClass[[i]], "list")
  expect_identical(p$default[[i]], list(p = 0.1))
})

test_that("the fit passes DEoptimControl, strategy and .c to runDEoptim", {
  ## The call needs a full simList, so this checks the parsed call instead of running it.
  exprs <- moduleSource()
  calls <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("runDEoptim"))) calls[[length(calls) + 1L]] <<- x
      lapply(as.list(x), walk)
    }
    invisible(NULL)
  }
  invisible(lapply(exprs, walk))
  expect_length(calls, 1L)
  args <- as.list(calls[[1L]])
  expect_identical(deparse(args$DEoptimControl), "P(sim)$DEoptimControl")
  expect_identical(deparse(args$strategy), "P(sim)$strategy")
  expect_identical(deparse(args$.c), "P(sim)$.c")
})
