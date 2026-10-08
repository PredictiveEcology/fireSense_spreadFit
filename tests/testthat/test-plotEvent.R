## The `plot` event (scheduled when `mode` has "visualize") draws one histogram per coefficient.
## It took the coefficient names from `sim$fireSense_SpreadFitted$bestCoef`, an object only the
## never-scheduled `retrieveDEOptim` event created. The names were therefore NULL, which is not an
## error: tidyr::gather() pooled every coefficient into a single unlabelled histogram.
## The names are `names(P(sim)$lower)`, which is what the `run` event uses for the same coefficients.

test_that("the plot event draws one histogram per coefficient, named as the run event names them", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("tidyr")
  rec <- new.env(); fitRec <- new.env()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE, SNLL_FS_thresh = 250L, mode = c("fit", "visualize")))
  mockFitAndLedger(sim, fitRec)
  mockInModule(sim, runSpreadWithoutDEoptim = function(...) NULL)
  sim <- suppressWarnings(suppressMessages(SpaDES.core::spades(sim)))
  expect_true("plot" %in% SpaDES.core::completed(sim)$eventType)
  lower <- SpaDES.core::params(sim)[[moduleName]]$lower
  expect_s3_class(sim$fsSpreadFit_hists, "ggplot")
  expect_setequal(unique(as.character(sim$fsSpreadFit_hists$data$key)), names(lower))
  expect_identical(nrow(sim$fsSpreadFit_hists$data), length(lower) * nrow(sim$DE[[1]]$member$pop))
})
