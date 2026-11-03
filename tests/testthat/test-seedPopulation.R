## DEoptim's own start is drawn blind: for an ELF whose feasible region is a sliver of the bounds it can have no
## member the objective accepts, and clusters stops the fit after generation 1. The fit then starts again from a
## gate-screened population (test-feasibleInitialPopulation.R). A population given as `initialpop` starts it first.

test_that("without an initialpop the first call has no starting population, and the retry has the screened one", {
  rec <- seededRun()
  expect_length(rec$calls, 1L)
  expect_null(rec$calls[[1]]$DEoptimControl$initialpop)
  rec <- seededRun(failFirst = TRUE)
  expect_length(rec$calls, 2L)
  expect_null(rec$calls[[1]]$DEoptimControl$initialpop)
  expect_equal(unname(rec$calls[[2]]$DEoptimControl$initialpop[1, ]), unname(rec$lower + 0.1))
})

test_that("a user-supplied initialpop reaches DEoptim first, then uniform draws inside the bounds", {
  ## runDEoptim() ignores its own initialpop argument: the given population never reached DEoptim (2026-10-26)
  first <- seededRun()
  given <- passingSets(first$lower, first$upper)[3:1, , drop = FALSE]
  rec <- seededRun(params = list(initialpop = given))
  pop <- rec$calls[[1]]$DEoptimControl$initialpop
  expect_false(is.null(pop))
  expect_equal(unname(pop[1:3, ]), unname(given))
  expect_identical(colnames(pop), names(first$lower))
  expect_gte(nrow(pop), 10L * length(first$lower))
  lo <- rep(first$lower, each = nrow(pop)); up <- rep(first$upper, each = nrow(pop))
  expect_true(all(pop >= lo & pop <= up))
})

test_that("an initialpop outside the bounds is an error, not a silent start outside them", {
  first <- seededRun()
  given <- passingSets(first$lower, first$upper)
  given[1, 1] <- first$upper[[1]] + 1
  expect_error(seededRun(params = list(initialpop = given)), "outside")
})

test_that("a cross-validation fold has no starting population", {
  rec <- new.env(); rec$calls <- list()
  sim <- toySim(list(heldOutFold = 1L, simulateMembers = 3L))
  mockFitAndLedger(sim, rec)
  mockInModule(sim,
    runDEoptim = function(...) { a <- list(...); rec$calls <- c(rec$calls, list(a)); toyDE(length(a$lower)) },
    Cache = function(FUN, pop, fnArgs, ...) {
      if (!is.function(FUN)) return(FUN)
      data.table::data.table(member = 1L, yr = names(fnArgs$historicalFires), rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
  suppressMessages(SpaDES.core::spades(sim))
  expect_null(rec$calls[[1]]$DEoptimControl$initialpop)
})

## 2026-10-26: a population given as the initialpop parameter (an emulator's best points) never reached DEoptim:
## the fit started from DEoptim's own random draw (generation 1 median at the fail value).
test_that("startPopulation() gives DEoptim the initialpop parameter when it is set", {
  lo <- c(a = 0, b = -1); up <- c(a = 1, b = 1)
  given <- matrix(c(0.1, 0.2, 0.3, 0, 0.5, -0.5), ncol = 2, dimnames = list(NULL, c("a", "b")))
  s <- startPopulation(given, 10, lo, up, 1)
  expect_equal(nrow(s), 10); expect_equal(unname(s[1:3, ]), unname(given))   # its rows first, then uniform draws
  expect_true(all(s[, "a"] >= 0 & s[, "a"] <= 1))
  expect_equal(unname(startPopulation(given, 2, lo, up, 1)), unname(given))  # n below its rows: all kept
  expect_identical(colnames(s), c("a", "b"))
  ## and reaches DEoptim's control
  expect_equal(deoptimControlStart(list(NP = 3), s)$initialpop, s)
  expect_error(startPopulation(given[, 2:1], 10, lo, up, 1), "not names\\(lower\\)")
  expect_error(startPopulation(given[, 1, drop = FALSE], 10, lo, up, 1), "columns")
  given[1, 1] <- 2
  expect_error(startPopulation(given, 10, lo, up, 1), "outside lower/upper")
  ## without it, NULL (DEoptim draws)
  expect_null(startPopulation(NULL, 10, lo, up, 1))
})
