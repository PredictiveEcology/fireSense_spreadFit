## The threshold calibration finds parameter sets that pass the threshold (ELF 6.6.1, weighted = "log": 6 of
## ~96 gate-passing sets). DEoptim's own start was drawn blind, twice had no member under the threshold, and
## clusters stopped the fit after generation 1. The fit now starts from the calibration's passing sets.

## three sets inside the bounds, columns as DEoptim uses them (the bounds are made by spreadFitPrepare, so the
## mocked calibration builds them from the bounds it is given)
passingSets <- function(lo, up) {
  m <- rbind(lo + 0.25 * (up - lo), lo + 0.5 * (up - lo), lo + 0.75 * (up - lo))
  colnames(m) <- names(lo)
  m
}

noFeasible <- function() structure(class = c("deoptimNoFeasibleMember", "error", "condition"),
                                   list(message = "every member of the population holds the fail value", call = NULL))

## The run event with the calibration and runDEoptim() replaced. `calibration` is what the cached
## calibration returns; `failFirst` makes the first fit stop as clusters does.
seededRun <- function(calibration, failFirst = FALSE, params = list()) {
  ## the module's Cache() of the calibration would serve the previous test's mocked result
  withr::local_options(reproducible.useCache = FALSE)
  rec <- new.env(); rec$calls <- list()
  sim <- toySim(c(list(firstBlockThresh = NULL), params))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runSpreadWithoutDEoptim = function(lower, upper, ...) {
      rec$lower <- lower; rec$upper <- upper; rec$passing <- passingSets(lower, upper)
      if (is.function(calibration)) calibration(rec) else calibration
    },
    runDEoptim = function(...) {
      args <- list(...)
      rec$calls <- c(rec$calls, list(args))
      if (failFirst && length(rec$calls) == 1L) stop(noFeasible())
      toyDE(length(args$lower))
    },
    feasibleInitialPopulation = function(n, lower, upper, gates, seed, ...) {
      m <- matrix(rep(lower, each = n) + 0.1, n, length(lower), dimnames = list(NULL, names(lower)))
      attr(m, "feasible") <- 5L
      m
    })
  suppressMessages(SpaDES.core::spades(sim))
  rec
}

test_that("the first fit starts from the calibration's passing sets, then uniform draws inside the bounds", {
  rec <- seededRun(function(r) list(thresh = 777, passingPars = r$passing))
  pop <- rec$calls[[1]]$DEoptimControl$initialpop
  expect_identical(unname(pop[1:3, ]), unname(rec$passing))
  expect_identical(colnames(pop), names(rec$lower))
  expect_gte(nrow(pop), 10L * length(rec$lower))
  lo <- rep(rec$lower, each = nrow(pop)); up <- rep(rec$upper, each = nrow(pop))
  expect_true(all(pop >= lo & pop <= up))
  ## the threshold is the number it always was
  expect_identical(as.vector(rec$calls[[1]]$thresh), 777)
  expect_null(attributes(rec$calls[[1]]$thresh)$passingPars)
})

test_that("the retry puts the passing sets first, before the gate-screened population", {
  rec <- seededRun(function(r) list(thresh = 777, passingPars = r$passing), failFirst = TRUE)
  expect_length(rec$calls, 2L)
  pop <- rec$calls[[2]]$DEoptimControl$initialpop
  expect_identical(unname(pop[1:3, ]), unname(rec$passing))
  expect_equal(unname(pop[4, ]), unname(rec$lower + 0.1))    # then the screened sets, as before
  expect_identical(nrow(pop), nrow(rec$calls[[1]]$DEoptimControl$initialpop))
})

test_that("a passing set outside the bounds is an error, not a silent start outside them", {
  expect_error(seededRun(function(r) { m <- r$passing; m[1, 1] <- r$upper[[1]] + 1; list(thresh = 1, passingPars = m) }),
               "outside")
})

test_that("with no passing sets the fit is as before: no starting population on the first call", {
  for (cal in list(list(thresh = 777, passingPars = NULL), 777)) {     # the second is a calibration cached before this change
    rec <- seededRun(cal)
    expect_null(rec$calls[[1]]$DEoptimControl$initialpop)
    expect_identical(as.vector(rec$calls[[1]]$thresh), 777)
  }
  rec <- seededRun(777, failFirst = TRUE)
  expect_equal(unname(rec$calls[[2]]$DEoptimControl$initialpop[1, ]), unname(rec$lower + 0.1))
})

test_that("a user-supplied initialpop reaches DEoptim, first, instead of the calibration's passing sets", {
  ## the old test expected no initialpop in DEoptimControl, relying on runDEoptim()'s own initialpop argument,
  ## which runDEoptim() ignores: the given population never reached DEoptim (2026-10-26)
  first <- seededRun(function(r) list(thresh = 777, passingPars = r$passing))
  given <- passingSets(first$lower, first$upper)[3:1, , drop = FALSE]
  rec <- seededRun(function(r) list(thresh = 777, passingPars = r$passing), params = list(initialpop = given))
  pop <- rec$calls[[1]]$DEoptimControl$initialpop
  expect_false(is.null(pop))
  expect_equal(unname(pop[1:3, ]), unname(given))
  expect_identical(colnames(pop), names(first$lower))
})

test_that("a cross-validation fold uses the full data's threshold and no starting population", {
  withr::local_options(reproducible.useCache = FALSE)
  rec <- new.env(); rec$calls <- list()
  sim <- toySim(list(heldOutFold = 1L, firstBlockThresh = NULL, simulateMembers = 3L))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runSpreadWithoutDEoptim = function(lower, upper, ...) list(thresh = 123, passingPars = passingSets(lower, upper)),
    runDEoptim = function(...) { a <- list(...); rec$calls <- c(rec$calls, list(a)); toyDE(length(a$lower)) },
    Cache = function(FUN, pop, fnArgs, ...) {
      if (!is.function(FUN)) return(FUN)
      data.table::data.table(member = 1L, yr = names(fnArgs$historicalFires), rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
  suppressMessages(SpaDES.core::spades(sim))
  expect_identical(as.vector(rec$calls[[1]]$thresh), 123)
  expect_null(rec$calls[[1]]$DEoptimControl$initialpop)
})

test_that("withPassingFirst() keeps the population's size and bounds, and is the identity without passing sets", {
  lo <- c(a = 0, b = 10); up <- c(a = 1, b = 20)
  pop <- matrix(c(0.1, 0.2, 0.3, 11, 12, 13), 3, dimnames = list(NULL, c("a", "b")))
  attr(pop, "feasible") <- 2L
  ## columns in another order than `lower`
  passing <- matrix(c(15, 19, 0.9, 0.8), 2, dimnames = list(NULL, c("b", "a")))
  out <- withPassingFirst(passing, pop, lo, up)
  expect_identical(dim(out), dim(pop))
  expect_equal(as.vector(out), c(0.9, 0.8, 0.1, 15, 19, 11))   # passing sets, then the population's first row
  expect_identical(colnames(out), c("a", "b"))
  expect_identical(attr(out, "feasible"), 3L)           # cut to the population's 3 rows
  expect_identical(withPassingFirst(NULL, pop, lo, up), pop)
})

## 2026-10-26: a population given as the initialpop parameter (an emulator's best points) never reached DEoptim:
## the fit started from DEoptim's own random draw (generation 1 median at the fail value).
test_that("startPopulation() gives DEoptim the initialpop parameter when it is set", {
  lo <- c(a = 0, b = -1); up <- c(a = 1, b = 1)
  given <- matrix(c(0.1, 0.2, 0.3, 0, 0.5, -0.5), ncol = 2, dimnames = list(NULL, c("a", "b")))
  s <- startPopulation(given, NULL, 10, lo, up, 1)
  expect_equal(nrow(s), 10); expect_equal(unname(s[1:3, ]), unname(given))   # its rows first, then uniform draws
  expect_true(all(s[, "a"] >= 0 & s[, "a"] <= 1))
  expect_equal(unname(startPopulation(given, NULL, 2, lo, up, 1)), unname(given))  # n below its rows: all kept
  ## it wins over passing sets
  expect_equal(unname(startPopulation(given, given[3:1, , drop = FALSE], 10, lo, up, 1)[1:3, ]), unname(given))
  ## and reaches DEoptim's control
  expect_equal(deoptimControlStart(list(NP = 3), s)$initialpop, s)
  expect_error(startPopulation(given[, 2:1], NULL, 10, lo, up, 1), "not names\\(lower\\)")
  ## without it: passing sets first, then uniform draws; without either, NULL (DEoptim draws)
  s <- startPopulation(NULL, given[1, , drop = FALSE], 10, lo, up, 1)
  expect_equal(nrow(s), 10); expect_equal(unname(s[1, ]), unname(given[1, ]))
  expect_null(startPopulation(NULL, NULL, 10, lo, up, 1))
})
