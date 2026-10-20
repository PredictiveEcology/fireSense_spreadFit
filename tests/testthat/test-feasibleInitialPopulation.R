## A fit starts from DEoptim's uniform draw between the (deliberately wide, sign-only) bounds. For an ELF
## whose feasible region is a sliver of that box, no member of the start passes the objective's
## spreadProb gates, every value is the fail value, and DE never improves (ELF 6.6.2, 2026-10-20:
## 0 of 2000 uniform draws passed; its fit died at generation 10). clusters now stops such a fit after
## its first generation with a `deoptimNoFeasibleMember` error; the module then starts again from a
## population drawn as the threshold calibration draws its trials (gate-screened), filled up with
## uniform draws. The first fit's cache key is unchanged, so fits in progress resume as before.

lower <- c(maxAsymptote = 0.14, CMD = 0, CMDsm = 0, youngAge = -50)
upper <- c(maxAsymptote = 0.27, CMD = 50, CMDsm = 50, youngAge = 0)

## gates that accept a thin slice of the box: both climate coefficients under 5
sliverGates <- function(parList, pass = TRUE) {
  ok <- vapply(parList, function(p) p[["CMD"]] < 5 && p[["CMDsm"]] < 5, logical(1))
  if (pass) ok else data.frame(pass = ok, burny = FALSE, notSpread = !ok, medianOK = TRUE)
}
noneGates <- function(parList, pass = TRUE) {
  if (pass) rep(FALSE, length(parList))
  else data.frame(pass = FALSE, burny = FALSE, notSpread = TRUE, medianOK = rep(c(TRUE, FALSE), length.out = length(parList)))
}

test_that("the population has n rows, the accepted sets first, all inside the bounds, and is reproducible", {
  expect_message(pop <- feasibleInitialPopulation(20L, lower, upper, gates = sliverGates, seed = 7L),
                 "Starting population: [0-9]+ of 20 members pass")
  expect_identical(dim(pop), c(20L, 4L))
  expect_identical(colnames(pop), names(lower))
  k <- attr(pop, "feasible")
  expect_gte(k, 4L)
  expect_true(all(sliverGates(lapply(seq_len(k), function(i) pop[i, ]))))
  expect_true(all(pop >= rep(lower, each = 20) & pop <= rep(upper, each = 20)))
  again <- suppressMessages(feasibleInitialPopulation(20L, lower, upper, gates = sliverGates, seed = 7L))
  expect_identical(unclass(again), unclass(pop))
  other <- suppressMessages(feasibleInitialPopulation(20L, lower, upper, gates = sliverGates, seed = 8L))
  expect_false(identical(unclass(other), unclass(pop)))
})

test_that("fewer than minFeasible accepted sets is an error that reports the gate counts", {
  expect_error(suppressMessages(feasibleInitialPopulation(8L, lower, upper, gates = noneGates, seed = 1L)),
               "only 0 of 800 parameter sets .* pass .* not spread out enough 200, median in range 100")
})

test_that("DEoptimControl is unchanged without a starting population, and carries one when given", {
  ctrl <- list(p = 0.1)
  expect_identical(deoptimControlStart(ctrl, NULL), ctrl)
  pop <- matrix(1, 2, 4)
  expect_identical(deoptimControlStart(ctrl, pop), list(p = 0.1, initialpop = pop))
  expect_identical(deoptimControlStart(NULL, pop), list(initialpop = pop))
})

test_that("spreadGates() applies fireSenseUtils::spreadProbGates() to the fit's data and settings", {
  seen <- new.env()
  testthat::local_mocked_bindings(
    spreadProbGates = function(par, ...) { seen$args <- list(...); data.frame(pass = c(TRUE, FALSE), burny = FALSE, notSpread = c(FALSE, TRUE), medianOK = TRUE) },
    .package = "fireSenseUtils")
  covs <- list(annualDTx1000 = "a", nonAnnualDTx1000 = "na", historicalFires = "hf")
  g <- spreadGates(covs, formulaToFit = "~ CMD", covMinMax = "mm", covCentre = NULL, climateRef = NULL,
                   mutuallyExclusive = list(youngAge = "x"), maxFireSpread = 0.276, escapeSizeHa = 50,
                   landscape = "ras", link = "powerLaw", powerLawGamma = 1.7)
  expect_identical(g(list(lower, upper)), c(TRUE, FALSE))
  expect_identical(seen$args$annualDTx1000, "a")
  expect_identical(seen$args$historicalFires, "hf")
  expect_identical(seen$args$link, "powerLaw")
  expect_identical(seen$args$powerLawGamma, 1.7)
  expect_false(seen$args$doAssertions)
  expect_identical(g(list(lower, upper), pass = FALSE)$notSpread, c(FALSE, TRUE))
})

noFeasible <- function() structure(class = c("deoptimNoFeasibleMember", "error", "condition"),
                                   list(message = "every member of the population holds the fail value", call = NULL))

## The run event with runDEoptim() replaced: the first call stops as clusters does, the second gets the
## gate-screened start in DEoptimControl and returns a fit.
retriedFit <- function(secondCall = function(args) toyDE(length(args$lower))) {
  rec <- new.env(); rec$calls <- list()
  sim <- toySim()
  mockFitAndLedger(sim, rec)
  start <- matrix(0.5, 90, 8); attr(start, "feasible") <- 5L
  mockInModule(sim,
    runDEoptim = function(...) {
      args <- list(...)
      rec$calls <- c(rec$calls, list(args))
      if (length(rec$calls) == 1L) stop(noFeasible())
      secondCall(args)
    },
    feasibleInitialPopulation = function(n, lower, upper, gates, seed, ...) {
      rec$n <- n; rec$seed <- seed; rec$gates <- gates
      start
    })
  rec$start <- start
  rec$sim <- sim
  rec
}

test_that("a fit clusters cannot start is run again from a gate-screened population, with the same other settings", {
  rec <- retriedFit()
  expect_message(sim <- SpaDES.core::spades(rec$sim), "Drawing a starting population that passes the spreadProb gates")
  expect_length(rec$calls, 2L)
  first <- rec$calls[[1]]; second <- rec$calls[[2]]
  ## the first call is the fit as always run: no starting population, so its cache key is unchanged
  expect_identical(first$DEoptimControl, list(p = 0.1))
  expect_null(first$DEoptimControl$initialpop)
  ## the retry carries the population and nothing else differs
  expect_identical(second$DEoptimControl, list(p = 0.1, initialpop = rec$start))
  expect_identical(second[setdiff(names(second), "DEoptimControl")], first[setdiff(names(first), "DEoptimControl")])
  ## drawn for at least 10 x parameters members, from the ELF's seed, with the fit's gates
  expect_identical(rec$n, 10L * length(first$lower))
  expect_identical(rec$seed, .elfSeed("9.9"))
  expect_true(is.function(rec$gates))
  expect_true(!is.null(sim$DE))
})

test_that("a second stop is final and names the population it started from", {
  rec <- retriedFit(secondCall = function(args) stop(noFeasible()))
  expect_error(suppressMessages(SpaDES.core::spades(rec$sim)),
               "failed again from a starting population of which 5 of 90 members pass")
  expect_length(rec$calls, 2L)
})
