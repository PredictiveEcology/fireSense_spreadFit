## The `debug` event called runSpreadWithoutDEoptim() with tests = "", which switches off
## every test in the objective function: no MAD, no adTest, no CRPS, so a debug run
## exercised none of the code it exists to debug. Every call must pass the `DEoptimTests`
## parameter, as the DEoptim fit does.
##
## The calls need a full simList, so this checks the module's source instead of running
## it: parsed, not sourced, so no module machinery is needed.
test_that("every runSpreadWithoutDEoptim() call passes the DEoptimTests parameter", {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_spreadFit.R"), keep.source = FALSE)

  calls <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("runSpreadWithoutDEoptim")))
        calls[[length(calls) + 1L]] <<- x
      for (i in seq_along(x)[-1L]) if (is.call(x[[i]])) walk(x[[i]])
    }
  }
  for (e in exprs) walk(e)

  expect_gte(length(calls), 1L) # the debug event
  for (cl in calls)
    expect_identical(deparse(cl$tests), "P(sim)$DEoptimTests")
})
