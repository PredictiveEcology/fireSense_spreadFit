## initialPopMethod = "emulator": a Gaussian-process emulator of the objective chooses DEoptim's starting population
## (prototype 2026-10-25/26: ELF 6.3.1 seeded reached in 16 generations the median the unseeded fit took 249 for).

test_that("with initialPopMethod = 'emulator' the emulator's points start DEoptim", {
  rec <- seededRun(params = list(initialPopMethod = "emulator"), mocks = emulatorMock)
  expect_false(is.null(rec$emu))
  call <- rec$calls[[1]]
  pop <- call$DEoptimControl$initialpop
  expect_identical(unname(pop[seq_len(nrow(rec$emuPop)), ]), unname(rec$emuPop))
  expect_identical(colnames(pop), names(rec$lower))
  ## the stopping rule is the one an unseeded fit gets
  expect_identical(call$convergence, clusters::deoptimConvergence(NULL))
  ## the population size and the cores are the fit's
  expect_identical(rec$emu$NP, fireSenseUtils::deoptimPopulationSize(length(rec$lower)))
  expect_identical(rec$emu$cores, call$cores)
  ## the objective the emulator evaluates is the fit's: the same arguments
  oa <- rec$emu$objArgs
  shared <- setdiff(intersect(names(oa), names(call)), "fireBufferedListDT")
  expect_true(all(c("formulaToFit", "covMinMax", "tests", "Nreps", "landscape", "link") %in% shared))
  for (nm in shared) expect_identical(oa[[nm]], call[[nm]], label = nm)
  expect_false("thresh" %in% names(oa))
  expect_false(any(c("itermax", "cores", "DEoptimControl", "lower") %in% names(oa)))
  expect_identical(oa$doAssertions, call$doObjFunAssertions)
  expect_identical(oa$fitYearSpreadSD, "yearSpreadSD" %in% names(rec$lower))
  expect_identical(oa$fireBufferedListDT, fireSenseUtils::addBufferEdge(call$fireBufferedListDT, call$landscape))
})

test_that("an emulator-started fit stops by the run's rule, as every fit does", {
  withr::local_options(clusters.deoptimConvergenceWindow = 300L, clusters.deoptimMinGenerations = 350L)
  rec <- seededRun(params = list(initialPopMethod = "emulator"), mocks = emulatorMock)
  expect_false(is.null(rec$emu))
  expect_identical(rec$calls[[1]]$convergence$window, 300L)
  expect_identical(rec$calls[[1]]$convergence$minGenerations, 350L)
})

test_that("a convergence parameter wins over the run's options for an emulator-started fit", {
  withr::local_options(clusters.deoptimConvergenceWindow = 300L, clusters.deoptimMinGenerations = 350L)
  rec <- seededRun(params = list(initialPopMethod = "emulator", convergence = list(window = 250L)),
                   mocks = emulatorMock)
  expect_identical(rec$calls[[1]]$convergence$window, 250L)
  expect_identical(rec$calls[[1]]$convergence$minGenerations, 350L)
})

test_that("initialPopMethod = 'calibration' (the default) is the fit as before: no emulator, the usual stopping rule", {
  rec <- seededRun(mocks = emulatorMock)
  expect_null(rec$emu)
  expect_identical(rec$calls[[1]]$convergence, clusters::deoptimConvergence(NULL))
  expect_null(rec$calls[[1]]$DEoptimControl$initialpop)
})

test_that("more parameters than emulatorMaxParams, or a given initialpop, start the fit as 'calibration' does", {
  rec <- seededRun(params = list(initialPopMethod = "emulator", emulatorMaxParams = 2L), mocks = emulatorMock)
  expect_null(rec$emu)
  expect_identical(rec$calls[[1]]$convergence, clusters::deoptimConvergence(NULL))
  expect_null(rec$calls[[1]]$DEoptimControl$initialpop)
  expect_message(expect_false(useEmulator(list(initialPopMethod = "emulator", emulatorMaxParams = 2L), 3L)),
                 "more than emulatorMaxParams")
  expect_true(useEmulator(list(initialPopMethod = "emulator", emulatorMaxParams = 18L), 3L))
  expect_false(useEmulator(list(initialPopMethod = "calibration", emulatorMaxParams = 18L), 3L))
  expect_error(useEmulator(list(initialPopMethod = "GP", emulatorMaxParams = 18L), 3L), "initialPopMethod")

  first <- seededRun()
  given <- passingSets(first$lower, first$upper)[3:1, , drop = FALSE]
  rec <- seededRun(params = list(initialPopMethod = "emulator", initialpop = given), mocks = emulatorMock)
  expect_null(rec$emu)
  expect_equal(unname(rec$calls[[1]]$DEoptimControl$initialpop[1:3, ]), unname(given))
})

test_that("the search space is the full bounds, signed log1p for wide one-sided parameters, and maps back exactly", {
  lo <- c(maxAsymptote = 0.14, `(Intercept)` = -50, CMD = 0, youngAge = -50, fuel = 0, nf = -50, sd = 0)
  up <- c(maxAsymptote = 0.265, `(Intercept)` = 50, CMD = 50, youngAge = 0, fuel = 100, nf = 50, sd = 10)
  sp <- emulatorSpace(lo, up, logSpanMin = 20)
  ## one bound 0 and the other at least 20 from it: log; a bound off 0, or a narrow one-sided range: linear
  expect_identical(sp$logPars, c("CMD", "youngAge", "fuel"))
  ## the bounds are the unit cube's faces
  expect_equal(unname(sp$toUnit(rbind(lo, up))), rbind(rep(0, 7), rep(1, 7)))
  P <- withr::with_seed(1, sapply(seq_along(lo), function(j) runif(50, lo[j], up[j])))
  colnames(P) <- names(lo)
  U <- sp$toUnit(P)
  expect_true(all(U >= 0 & U <= 1))
  expect_equal(sp$fromUnit(U), P)                                          # round trip
  expect_identical(colnames(sp$fromUnit(U)), names(lo))
  back <- sp$fromUnit(rbind(rep(0, 7), rep(1, 7)))
  expect_true(all(t(back) >= lo & t(back) <= up))                         # never outside the bounds
  ## half the unit width of a log parameter is near its 0 bound: fuel 0-9 and youngAge -6.1-0, not 0-50 and -25-0
  mid <- sp$fromUnit(matrix(0.5, 1, 7))
  expect_equal(mid[1, c("fuel", "youngAge", "CMD")], c(fuel = expm1(log1p(100) / 2), youngAge = -expm1(log1p(50) / 2),
                                                       CMD = expm1(log1p(50) / 2)))
  expect_equal(mid[1, c("(Intercept)", "nf", "sd")], c(`(Intercept)` = 0, nf = 0, sd = 5))
  ## logSpanMin = Inf: all linear
  expect_identical(emulatorSpace(lo, up, logSpanMin = Inf)$logPars, character(0))
})

test_that("the trust region's first width is the best points' box, widened; the box is centred on the best point", {
  X <- rbind(c(0.40, 0.10), c(0.50, 0.30), c(0.90, 0.90), c(0.10, 0.95), c(0.45, 0.20),
             c(0.20, 0.80), c(0.70, 0.60), c(0.30, 0.50), c(0.60, 0.40), c(0.80, 0.70))
  score <- c(1, 2, 9, 8, 3, 7, 6, 5, 4, 10)
  noMargin <- function(b) { attr(b, "margin") <- NULL; b }
  w <- emulatorTrustWidth(X, score, topFrac = 0.3, widen = 0.5)          # the 3 best: rows 1, 2, 5
  expect_equal(w, c(0.2, 0.4))                                           # x 0.4-0.5, y 0.1-0.3, each doubled
  expect_equal(emulatorTrustWidth(X, score, topFrac = 0.3, widen = 0), c(0.1, 0.2))
  box <- emulatorTrustBox(X[1, ], w, widen = 0.5)
  expect_equal(noMargin(box), rbind(lower = c(0.3, 0), upper = c(0.5, 0.3)))   # y clipped at 0
  expect_equal(attr(box, "margin"), c(0.05, 0.1))                        # what widen added, before clipping
  expect_equal(noMargin(emulatorTrustBox(X[1, ], w / 2, widen = 0.5)), rbind(lower = c(0.35, 0), upper = c(0.45, 0.2)))
  expect_equal(noMargin(emulatorTrustBox(c(0.5, 0.5), c(4, 0.2), widen = 0.5)), rbind(lower = c(0, 0.4), upper = c(1, 0.6)))
})

test_that("a trust round that improves by more than the GP's uncertainty shrinks the box; repeated failures stop the search", {
  ctl <- emulatorControlOf(list(trustShrink = 0.5, trustExpand = 1.5, trustSuccessTol = 1L, trustFailureTol = 2L))
  st <- emulatorTrustStep(list(scale = 1, successes = 0L, failures = 0L), improvement = 30, tol = 9, ctl)
  expect_equal(st$scale, 0.5); expect_false(st$stop)
  st <- emulatorTrustStep(st, improvement = 12, tol = 9, ctl)
  expect_equal(st$scale, 0.25)
  st <- emulatorTrustStep(st, improvement = 3, tol = 9, ctl)              # a failure: the box grows a little
  expect_equal(st$scale, 0.375); expect_identical(st$failures, 1L); expect_false(st$stop)
  st <- emulatorTrustStep(st, improvement = 20, tol = 9, ctl)             # a success resets the failures
  expect_identical(st$failures, 0L); expect_equal(st$scale, 0.1875)
  st <- emulatorTrustStep(st, improvement = 1, tol = 9, ctl)
  st <- emulatorTrustStep(st, improvement = 1, tol = 9, ctl)              # two failures in a row: stop
  expect_true(st$stop); expect_identical(st$reason, "improvement below noise")
  ## the box keeps its own scale, which can grow past its first width (emulatorTrustBox() keeps it in the cube)
  st <- emulatorTrustStep(list(scale = 0.9, successes = 0L, failures = 0L), improvement = 0, tol = 9, ctl)
  expect_equal(st$scale, 1.35)
  ## two successes needed before it shrinks
  ctl2 <- emulatorControlOf(list(trustSuccessTol = 2L))
  st <- emulatorTrustStep(list(scale = 1, successes = 0L, failures = 0L), improvement = 30, tol = 9, ctl2)
  expect_equal(st$scale, 1)
  expect_equal(emulatorTrustStep(st, improvement = 30, tol = 9, ctl2)$scale, ctl2$trustShrink)
})

test_that("a trust round whose best point is at the box's edge grows the box", {
  ## before, every success halved the box, so a search still moving toward a minimum outside the box shrank
  ## it below trustMinWidth, or to steps smaller than the noise, and stopped 7-15 noise SDs above the minimum
  box <- structure(rbind(lower = c(0.2, 0), upper = c(0.6, 0.4)), margin = c(0.1, 0.1))
  expect_false(emulatorAtBoxEdge(c(0.4, 0.2), box))
  expect_true(emulatorAtBoxEdge(c(0.25, 0.2), box))                      # in the lower margin of x
  expect_true(emulatorAtBoxEdge(c(0.4, 0.35), box))                      # in the upper margin of y
  expect_true(emulatorAtBoxEdge(c(0.7, 0.2), box))                       # outside the box
  expect_false(emulatorAtBoxEdge(c(0.4, 0.05), box))                     # y's lower side is the cube's face
  ## an upper side at the cube's face is not an edge either
  boxUp <- structure(rbind(lower = c(0.2, 0.5), upper = c(0.6, 1)), margin = c(0.1, 0.1))
  expect_false(emulatorAtBoxEdge(c(0.4, 0.95), boxUp))
  expect_true(emulatorAtBoxEdge(c(0.4, 0.55), boxUp))
  ## a box clipped on one side keeps the margin it was given (not its clipped width * widen / (1 + 2 * widen))
  clipped <- emulatorTrustBox(c(0.45, 0.1), c(0.2, 0.3), widen = 0.5)    # y: -0.05-0.25, lower side clipped
  expect_equal(attr(clipped, "margin"), c(0.05, 0.075))
  expect_equal(clipped[, 2L], c(lower = 0, upper = 0.25))
  expect_true(emulatorAtBoxEdge(c(0.45, 0.18), clipped))                 # in the 0.075 margin of the open side
  expect_false(emulatorAtBoxEdge(c(0.45, 0.17), clipped))
  ctl <- emulatorControlOf(list(trustShrink = 0.5, trustExpand = 1.25, trustSuccessTol = 1L))
  st <- emulatorTrustStep(list(scale = 0.5, successes = 0L, failures = 1L), improvement = 30, tol = 9, ctl,
                          atEdge = TRUE)
  expect_equal(st$scale, 0.625); expect_identical(st$failures, 0L); expect_false(st$stop)
  expect_equal(emulatorTrustStep(st, improvement = 30, tol = 9, ctl, atEdge = FALSE)$scale, 0.3125)
})

test_that("a search moving toward a minimum outside its first box grows the box while it is at the edge", {
  ## no GP and no noise: each round's best point is the true minimum clipped to the round's box, centred on the
  ## best point so far
  f <- function(P) rowSums(sweep(P, 2L, c(0.9, 0.85))^2)
  ctl <- emulatorControlOf(list(trustShrink = 0.5, trustExpand = 2, trustSuccessTol = 1L))
  X <- rbind(c(0.10, 0.10), c(0.15, 0.20), c(0.20, 0.15), c(0.25, 0.25), c(0.30, 0.10), c(0.10, 0.30))
  width <- emulatorTrustWidth(X, f(X), ctl$localTopFrac, ctl$localWiden)
  st <- list(scale = 1, successes = 0L, failures = 0L)
  scales <- numeric(0); edges <- logical(0)
  for (i in 1:4) {
    box <- emulatorTrustBox(X[which.min(f(X)), ], st$scale * width, ctl$localWiden)
    best <- pmin(pmax(c(0.9, 0.85), box[1L, ]), box[2L, ])
    edges <- c(edges, emulatorAtBoxEdge(best, box))
    before <- min(f(X)); X <- rbind(X, best)
    st <- emulatorTrustStep(st, before - min(f(X)), tol = 1e-6, ctl, utils::tail(edges, 1L))
    scales <- c(scales, st$scale)
  }
  expect_true(all(edges))
  expect_equal(scales, c(2, 4, 8, 16))                                   # the first width is 0.1 by 0.3
  expect_equal(min(f(X)), 0)                                             # the fourth box reaches the minimum
})

test_that("the starting population is good evaluated points, picked best first at least a minimum distance apart", {
  ## two good regions: a tight cluster of the 5 best around (0.2, 0.2), and the next 3 around (0.8, 0.8)
  X <- rbind(c(0.20, 0.20), c(0.21, 0.20), c(0.20, 0.21), c(0.22, 0.22), c(0.19, 0.19),
             c(0.80, 0.80), c(0.81, 0.80), c(0.80, 0.81), c(0.50, 0.10), c(0.50, 0.90))
  value <- c(1:8, 50, 60)
  pick <- emulatorSpreadPopulation(X, value, theta = c(0.1, 0.1), NP = 3L, poolFrac = 0.8, minDist = 1)
  ## the best, then the best of the other region; (0.21, 0.20) is 0.1 lengthscales from the best, so skipped
  expect_identical(pick[1:2], c(1L, 6L))
  expect_length(pick, 3L)
  expect_false(anyDuplicated(pick) > 0)
  expect_true(all(pick %in% 1:8))                                          # only from the pool, the best 80%
  ## no third point is 1 or 0.5 lengthscales from both; at 0.25, (0.22, 0.22), 0.28 from the best
  expect_identical(pick[3], 4L)
  ## NP rows even when the pool is smaller than NP: the pool grows to NP, then the distance relaxes to 0
  expect_identical(sort(emulatorSpreadPopulation(X, value, 0.1, NP = 10L, poolFrac = 0.1, minDist = 1)), 1:10)
  expect_length(emulatorSpreadPopulation(X[1:3, ], value[1:3], 0.1, NP = 8L, poolFrac = 0.2, minDist = 1), 3L)
})

test_that("emulatorControl names are checked", {
  expect_identical(emulatorControlOf(list(maxEvals = 10L))$maxEvals, 10L)
  expect_identical(emulatorControlOf()$designPerParam, emulatorDefaults()$designPerParam)
  expect_error(emulatorControlOf(list(maxEval = 10L)), "unknown element")
})

test_that("the GP is fitted to at most gpMaxPoints evaluations, always the best, and Inf keeps every one", {
  set.seed(3)
  U <- matrix(stats::runif(300 * 3), ncol = 3)
  U <- rbind(U, U[rep(1:2, each = 5), ])                                 # two points with 5 more replicates each
  y <- 100 + 50 * rowSums((U - 0.5)^2) + stats::rnorm(nrow(U))
  y[7] <- NA                                                             # a failed evaluation is the worst
  rows <- emulatorTrainRows(U, y, maxPoints = 80L, bestFrac = 0.5)
  expect_lte(length(rows), 80L)
  expect_gt(length(rows), 70L)                                           # the cap is filled, not just respected
  expect_identical(rows, sort(unique(rows)))
  ## the best distinct points by mean are all in, with their replicates
  keys <- emulatorPointKeys(U)
  m <- tapply(ifelse(is.na(y), Inf, y), keys, mean)
  bestKeys <- names(sort(m))[1:20]
  expect_true(all(which(keys %in% bestKeys) %in% rows))
  expect_false(7L %in% rows)
  ## a replicated point is in with all its rows or not at all
  for (k in unique(keys[301:310])) expect_true(all(which(keys == k) %in% rows) || !any(which(keys == k) %in% rows))
  ## the rest are spread out: the thinned points are farther apart than a random draw of the same size
  rest <- setdiff(rows, which(keys %in% names(sort(m))[1:40]))
  nn <- function(X) mean(apply(as.matrix(stats::dist(X)) + diag(Inf, nrow(X)), 1, min))
  expect_gt(nn(U[rest, ]), nn(U[sample(setdiff(seq_len(nrow(U)), rows), length(rest)), ]))
  ## Inf, or a cap above the evaluations, is every row: the earlier behaviour
  expect_identical(emulatorTrainRows(U, y, Inf, 0.5), seq_len(nrow(U)))
  expect_identical(emulatorTrainRows(U, y, nrow(U), 0.5), seq_len(nrow(U)))
  expect_identical(emulatorControlOf(list(gpMaxPoints = Inf))$gpMaxPoints, Inf)
})

## A noisy quadratic with its minimum near the 0 bound of two wide one-sided parameters, as fuel coefficients
## (bounds 0-100, fitted 0.5-2) and youngAge (-50-0, fitted -1) are; the objective changes by ~10 per unit of them
quadLower <- c(a = 0, b = -5, fuel = 0, young = -50); quadUpper <- c(a = 1, b = 5, fuel = 100, young = 0)
quadOpt <- c(a = 0.3, b = 1, fuel = 0.7, young = -1)
quadValue <- function(P)
  100 + 10 * rowSums(sweep(sweep(P[, names(quadLower), drop = FALSE], 2, quadOpt), 2, c(0.1, 1, 1, 1), "/")^2)
quadEmulate <- function(ctl, NP = 8L) {
  n <- 0L
  noisy <- function(P, seed) {
    expect_identical(colnames(P), names(quadLower))
    n <<- n + nrow(P)
    withr::with_seed(seed, quadValue(P) + stats::rnorm(nrow(P), 0, 2))
  }
  res <- suppressMessages(emulateObjective(noisy, quadLower, quadUpper, NP = NP, control = ctl, seed = 1L))
  res$nEvaluated <- n
  res
}
quadCtl <- list(designPerParam = 10L, globalMaxEvals = 150L, maxEvals = 400L, minEvals = 100L, batchSize = 12L,
                nCandidates = 500L, nPerturb = 300L, gpThreads = 1L)

test_that("the emulator fits its GP to the capped rows only, and to every row with gpMaxPoints = Inf", {
  fitted <- integer(0)
  testthat::local_mocked_bindings(fitEmulatorGP = function(U, z, ctl) { fitted <<- c(fitted, nrow(U)); stop("stop here") })
  for (cap in c(40L, Inf)) {
    fitted <- integer(0)
    expect_error(quadEmulate(utils::modifyList(quadCtl, list(gpMaxPoints = cap))), "stop here")
    expect_identical(fitted, if (is.finite(cap)) 40L else 64L)
  }
})

test_that("the emulator finds a minimum near the sign constraints of wide one-sided parameters within its budget", {
  ## the shape of a spread fit (ELF 6.6.1/6.6.2 bounds): the optimum's fuel coefficients, CMD and youngAge are 0.5-1.5
  ## against bounds of 50-100 (before, the search was linear in a box near the bounds and stopped ~80% above the
  ## optimum: this seed scored 144 here after 400 evaluations)
  lo <- c(maxA = 0.14, int = -50, CMD = 0, young = -50, dom = 0, sec = 0, nf = -50, sd = 0)
  up <- c(maxA = 0.265, int = 50, CMD = 50, young = 0, dom = 100, sec = 100, nf = 50, sd = 10)
  opt <- c(maxA = 0.25, int = 0.5, CMD = 1.4, young = -0.5, dom = 0.6, sec = 1.5, nf = -40, sd = 2.5)
  value <- function(P)
    100 + 10 * rowSums(sweep(sweep(P[, names(lo), drop = FALSE], 2, opt), 2, c(0.01, 1, 1, 1, 1, 1, 5, 1), "/")^2)
  n <- 0L
  noisy <- function(P, seed) {
    expect_identical(colnames(P), names(lo))
    n <<- n + nrow(P)
    withr::with_seed(seed, value(P) + stats::rnorm(nrow(P), 0, 2))
  }
  ctl <- list(designPerParam = 10L, globalMaxEvals = 250L, maxEvals = 400L, minEvals = 100L, batchSize = 12L,
              nCandidates = 1000L, nPerturb = 500L, gpThreads = 1L)
  res <- suppressMessages(emulateObjective(noisy, lo, up, NP = 8L, control = ctl, seed = 1L))
  expect_identical(colnames(res$pop), names(lo))
  expect_identical(nrow(res$pop), 8L)
  expect_false(anyDuplicated(res$pop) > 0)
  expect_true(all(t(res$pop) >= lo & t(res$pop) <= up))
  expect_lte(n, ctl$maxEvals + ctl$batchSize)
  expect_identical(nrow(res$evals), n)
  ## the population is evaluated points, with the mean of their evaluations
  key <- function(M) apply(signif(M[, names(lo), drop = FALSE], 10), 1, paste, collapse = ",")
  expect_true(all(key(res$pop) %in% key(as.matrix(res$evals[, names(lo), with = FALSE]))))
  ## the far corner is ~10^7 above the minimum, the design's best ~10^3
  expect_lt(min(value(res$pop)) - 100, 10)
})

test_that("the local phase is trust rounds that shrink around the minimum and stop when they improve less than the GP's uncertainty", {
  ## one global batch: the rounds find the minimum (over seeds 1..64 the gap is 0.38 at the median, >= 5 for 3,
  ## one of them a local GP flat at its cap; before the trust region kept its own scale, >= 5 for 19). This seed
  ## ends 0.39 from it with OpenBLAS, 0.22 with reference BLAS; the edge-growth tests above are the regression tests
  res <- quadEmulate(utils::modifyList(quadCtl, list(globalMaxEvals = 76L, maxEvals = 2000L)))
  st <- res$status[which(res$status$phase == "local"), ]
  expect_gte(max(st$round), 3L)
  expect_identical(res$stopReason, "improvement below noise")
  expect_lt(res$nEvaluated, 800L)                                      # 628 at most over seeds 1..64; 436 here
  expect_lt(min(st$scale), 0.5)
  expect_lt(min(quadValue(res$pop)) - 100, 5)
  ## two global batches find the minimum (63 of seeds 1..64 end within 5 of it): the search stops, near it
  res <- quadEmulate(utils::modifyList(quadCtl, list(globalMaxEvals = 88L, maxEvals = 2000L)))
  expect_identical(res$stopReason, "improvement below noise")
  expect_lt(res$nEvaluated, 600L)                                      # 472 at most over seeds 1..64; 136-232 here
  expect_lt(min(quadValue(res$pop)) - 100, 5)
})

test_that("the emulator's evaluations run on a cluster built by clusters::clusterSetup() from the fit's cores", {
  rec <- new.env()
  testthat::local_mocked_bindings(.package = "clusters", clusterSetup = function(...) {
    a <- list(...)
    rec$setup <- a
    rec$shipped <- mget(unlist(a$objsNeeded), envir = a$envir)
    list(cluster = list("w1", "w2", "w3"))
  })
  testthat::local_mocked_bindings(emulateObjective = function(evaluate, lower, upper, NP, control, seed) {
    rec$batchSize <- control$batchSize
    list(pop = NULL)
  })
  objArgs <- list(landscape = "L", annualDTx1000 = 1, nonAnnualDTx1000 = 2, fireBufferedListDT = 3,
                  historicalFires = 4, mutuallyExclusive = list(), tests = "adTest")
  lo <- c(a = 0); up <- c(a = 1)
  suppressMessages(emulatorStartPopulation(objArgs, lo, up, NP = 30L, control = list(), seed = 1L,
                                           cores = c("hostA", "hostB"), nCoresNeeded = 30L, runName = "toy"))
  expect_identical(rec$setup$cores, c("hostA", "hostB"))
  expect_identical(rec$setup$nCoresNeeded, 30L)
  expect_identical(rec$shipped, objArgs[emulatorDataArgs])
  expect_identical(rec$batchSize, 3L)                                   # one evaluation per worker per batch
  ## without cores, nothing is built and the evaluations run here
  rec$setup <- NULL
  suppressMessages(emulatorStartPopulation(objArgs, lo, up, NP = 30L, control = list(), seed = 1L, cores = NA))
  expect_null(rec$setup)
  expect_null(rec$batchSize)
})
