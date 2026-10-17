## The DEoptim early-stop settings are a fit setting the fit's cache key sees: the module resolves them once
## (parameter `convergence`, else the clusters options, else clusters' defaults) to the full list and hands it
## to runDEoptim() inside the Cache() call. A fit stopped under one rule, rerun under a stricter one, is a new
## call, which clusters resumes from its cached generations.

convergenceSeen <- function(params = list()) {
  rec <- new.env()
  sim <- toySim(c(list(doFit = NA), params))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runDEoptim = function(...) {
      a <- list(...)
      rec$convergence <- a$convergence
      toyDE(length(a$lower))
    },
    Cache = function(FUN, ..., omitArgs = NULL) {
      if (!is.function(FUN)) return(FUN)
      data.table::data.table(member = 1L, yr = "year2001", rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
  suppressMessages(SpaDES.core::spades(sim))
  rec$convergence
}

noConvergenceOptions <- function(env = parent.frame())
  withr::local_options(clusters.deoptimConvergenceWindow = NULL, clusters.deoptimMinGenerations = NULL,
                       clusters.deoptimConvergenceSE = NULL, clusters.deoptimConvergenceSDChange = NULL,
                       .local_envir = env)

test_that("convergence is NULL by default and the fit gets clusters' full default list", {
  noConvergenceOptions()
  expect_identical(convergenceSeen(), clusters::deoptimConvergence())
  expect_identical(convergenceSeen(), list(window = 300L, minGenerations = 350L, se = 1, sdChange = 0.1))
})

test_that("the parameter, then the clusters options, decide what the fit gets", {
  noConvergenceOptions()
  expect_identical(convergenceSeen(list(convergence = list(minGenerations = 800L))),
                   list(window = 300L, minGenerations = 800L, se = 1, sdChange = 0.1))
  withr::local_options(clusters.deoptimMinGenerations = 500L, clusters.deoptimConvergenceSE = 2)
  expect_identical(convergenceSeen(),
                   list(window = 300L, minGenerations = 500L, se = 2, sdChange = 0.1))
  expect_identical(convergenceSeen(list(convergence = list(minGenerations = 800L))),
                   list(window = 300L, minGenerations = 800L, se = 2, sdChange = 0.1))
})
