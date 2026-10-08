## After a fit, and for each held-out fold, the module compares the fit with its data: observed
## against simulated burned share by covariate, and the response curves, with Plots(). The
## simulation behind both is made once per fit, and only when `.plots` asks for a plot.
## R/validationFigures.R spreadFitValidationFigures().

## records the validation data calls and the Plots() calls instead of simulating and drawing
mockValidation <- function(sim, rec) {
  rec$data <- list(); rec$plots <- list()
  mockInModule(sim,
    spreadFitValidationData = function(par, seed, ...) {
      a <- list(...)
      d <- data.table::data.table(year = names(a$historicalFires)[1], pixelID = 1L, observed = 1L,
                                  simulated = 0.5, p = 0.2)
      rec$data[[length(rec$data) + 1L]] <- list(par = par, seed = seed, args = a, d = d)
      d
    },
    ## spreadFitPrep() also calls Plots() (covariate histograms); only these figures are recorded
    Plots = function(data, fn, filename, ...) {
      if (!startsWith(filename, "spreadFit")) return(invisible(NULL))
      a <- list(...)
      rec$plots[[length(rec$plots) + 1L]] <- list(data = data, fn = fn, types = a$types, path = a$path,
                                                  filename = filename)
      invisible(NULL)
    })
}

fitWithPlots <- function(params = list()) {
  rec <- new.env()
  sim <- toySim(c(list(stopIfNoPreRunFit = FALSE, profileReps = 0L, simulateMembers = 0L), params))
  mockFitAndLedger(sim, rec)
  mockValidation(sim, rec)
  list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
}

test_that("after a fit, .plots = 'png' makes both figures from one simulation of the best member", {
  out <- fitWithPlots(list(.plots = "png"))
  rec <- out$rec; sim <- out$sim
  expect_length(rec$data, 1L)                                   # data built once
  v <- rec$data[[1]]
  ## the best member by replicated mean: toyDE's member 7 has the lowest value (1)
  expect_identical(names(v$par), names(SpaDES.core::params(sim)[[moduleName]]$lower))
  expect_equal(unname(v$par[1]), 7)
  expect_identical(v$seed, .elfSeed("9.9"))
  ## the fit's simulation settings: every year, its replicates and its escape rule
  expect_identical(names(v$args$historicalFires), c("year2001", "year2002"))
  expect_identical(v$args$Nreps, 50L)
  expect_identical(v$args$escapeSizeHa, 50)
  expect_identical(v$args$jumpTries, 20)
  expect_identical(v$args$jumpMeanDist, 3)
  ## two Plots() calls on that one table, png, under figurePath(sim)
  expect_length(rec$plots, 2L)
  expect_identical(rec$plots[[1]]$data, v$d)
  expect_identical(rec$plots[[2]]$data, v$d)
  expect_identical(rec$plots[[1]]$fn, fireSenseUtils::plotSpreadFitValidation)
  expect_identical(rec$plots[[2]]$fn, fireSenseUtils::plotSpreadFitResponse)
  expect_identical(vapply(rec$plots, `[[`, "", "types"), c("png", "png"))
  ## figurePath(sim) inside the module's event: <outputPath>/figures/<module>
  expect_identical(vapply(rec$plots, `[[`, "", "path"),
                   rep(file.path(SpaDES.core::outputPath(sim), "figures", moduleName), 2))
  expect_identical(vapply(rec$plots, `[[`, "", "filename"),
                   c("spreadFitObservedVsSimulated_toyRun", "spreadFitResponseCurves_toyRun"))
})

test_that("with .plots NULL (the default) nothing is simulated or plotted", {
  expect_null(SpaDES.core::params(toySim())[[moduleName]]$.plots)
  rec <- fitWithPlots()$rec
  expect_length(rec$data, 0L)
  expect_length(rec$plots, 0L)
})

test_that("a held-out fold's figures score only its held-out years, from that fold's fit", {
  rec <- new.env()
  sim <- toySim(list(heldOutFold = 1L, simulateMembers = 3L, .plots = "png"))
  mockFitAndLedger(sim, rec)
  mockValidation(sim, rec)
  mockInModule(sim, Cache = function(FUN, pop, fnArgs, ...) {
    if (!is.function(FUN)) return(FUN)
    rec$heldOutArgs <- fnArgs
    data.table::data.table(member = 1L, yr = names(fnArgs$historicalFires), rep = 1L, ids = 1L, size = 4, sim = 8L)
  })
  suppressMessages(SpaDES.core::spades(sim))
  expect_length(rec$data, 1L)
  expect_identical(names(rec$data[[1]]$args$historicalFires), "year2001")   # fold 1 is year2001
  expect_identical(vapply(rec$plots, `[[`, "", "filename"),
                   c("spreadFitObservedVsSimulated_toyRun_cvFold1_heldOut",
                     "spreadFitResponseCurves_toyRun_cvFold1_heldOut"))
  ## the held-out simulation now keeps the fit's escape rule
  expect_identical(rec$heldOutArgs$escapeSizeHa, 50)
  expect_identical(rec$heldOutArgs$jumpTries, 20)
})
