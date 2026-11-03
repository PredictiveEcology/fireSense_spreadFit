## drawActivePars() draws from the logistic's active range: draws that saturate (fail the objective's
## spreadProb gates) are rejected, so what is returned are sets the objective will score.
## Here only a > 0.97 passes (3% of the uniform draws).

lower <- c(a = 0, b = 10)
upper <- c(a = 1, b = 20)
active <- function(par) par[["a"]] > 0.97
gates <- function(parList) vapply(parList, active, logical(1))

test_that("rejection sampling finds the sets that pass, named as the bounds", {
  set.seed(42L)
  d <- drawActivePars(4L, lower, upper, maxDraws = 400L, gates)
  expect_equal(d$accepted, 4L)
  expect_true(all(vapply(d$pars, active, logical(1))))
  expect_identical(names(d$pars[[1]]), names(lower))
  expect_gt(d$drawn, 4L)
})

test_that("the cap on draws is respected", {
  set.seed(42L)
  a40 <- vapply(1:40, function(i) runif(2, c(0, 10), c(1, 20))[1], numeric(1))
  cutoff <- sort(a40, decreasing = TRUE)[2]           # exactly two of the first 40 draws are >= this
  seen <- 0L
  two <- function(parList) { seen <<- seen + length(parList); vapply(parList, function(p) p[["a"]] >= cutoff, logical(1)) }
  set.seed(42L)
  d <- drawActivePars(4L, lower, upper, maxDraws = 40L, two)
  expect_equal(seen, 40L)                              # never gates more than the cap
  expect_equal(d$accepted, 2L)
  expect_equal(d$drawn, 40L)
  seen <- 0L
  set.seed(42L)
  d <- drawActivePars(4L, lower, upper, maxDraws = 40L, function(parList) { seen <<- seen + length(parList); rep(FALSE, length(parList)) })
  expect_equal(seen, 40L)
  expect_equal(d$accepted, 0L)
  expect_length(d$pars, 0L)
})

test_that("with a fixed seed the accepted draws are reproducible, whether the gates run in one process or two", {
  set.seed(7L); d1 <- drawActivePars(4L, lower, upper, 400L, gates, cores = 1L)
  set.seed(7L); d2 <- drawActivePars(4L, lower, upper, 400L, gates, cores = 2L)
  expect_identical(d2, d1)
  set.seed(8L); d3 <- drawActivePars(4L, lower, upper, 400L, gates, cores = 1L)
  expect_false(identical(d3$pars, d1$pars))
})
