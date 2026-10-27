## initialPopMethod = "emulator": a Gaussian-process emulator of the objective chooses DEoptim's starting population
## (prototype 2026-10-25/26: ELF 6.3.1 seeded reached in 16 generations the median the unseeded fit took 249 for).

## The emulator's search, replaced by a recorder that returns NP points at 0.3 of the way from lower to upper
emulatorMock <- function(rec) list(
  emulatorStartPopulation = function(objArgs, lower, upper, NP, control, seed, cores, nCoresNeeded, ...) {
    rec$emu <- list(objArgs = objArgs, NP = NP, cores = cores, nCoresNeeded = nCoresNeeded)
    pop <- matrix(rep(lower + 0.3 * (upper - lower), each = NP), NP, dimnames = list(NULL, names(lower)))
    rec$emuPop <- pop
    list(pop = pop)
  })
withPassing <- function(r) list(thresh = 777, passingPars = r$passing)

test_that("with initialPopMethod = 'emulator' the emulator's points start DEoptim, and the fit stops sooner", {
  rec <- seededRun(withPassing, params = list(initialPopMethod = "emulator"), mocks = emulatorMock)
  expect_false(is.null(rec$emu))
  call <- rec$calls[[1]]
  pop <- call$DEoptimControl$initialpop
  expect_identical(unname(pop[seq_len(nrow(rec$emuPop)), ]), unname(rec$emuPop))
  expect_identical(colnames(pop), names(rec$lower))
  ## only minGenerations and window change
  expect_identical(call$convergence, clusters::deoptimConvergence(list(window = 200L, minGenerations = 200L)))
  ## the population size and the cores are the fit's
  expect_identical(rec$emu$NP, fireSenseUtils::deoptimPopulationSize(length(rec$lower)))
  expect_identical(rec$emu$cores, call$cores)
  ## the objective the emulator evaluates is the fit's: the same arguments, with no early stop
  oa <- rec$emu$objArgs
  shared <- setdiff(intersect(names(oa), names(call)), c("fireBufferedListDT", "thresh"))
  expect_true(all(c("formulaToFit", "covMinMax", "tests", "Nreps", "landscape", "link") %in% shared))
  for (nm in shared) expect_identical(oa[[nm]], call[[nm]], label = nm)
  expect_identical(oa$thresh, Inf)
  expect_false(any(c("itermax", "cores", "DEoptimControl", "lower") %in% names(oa)))
  expect_identical(oa$doAssertions, call$doObjFunAssertions)
  expect_identical(oa$fitYearSpreadSD, "yearSpreadSD" %in% names(rec$lower))
  expect_identical(oa$fireBufferedListDT, fireSenseUtils::addBufferEdge(call$fireBufferedListDT, call$landscape))
})

test_that("a convergence parameter wins over the seeded fit's defaults", {
  rec <- seededRun(withPassing, params = list(initialPopMethod = "emulator", convergence = list(window = 250L)),
                   mocks = emulatorMock)
  expect_identical(rec$calls[[1]]$convergence, clusters::deoptimConvergence(list(window = 250L, minGenerations = 200L)))
})

test_that("initialPopMethod = 'calibration' (the default) is the fit as before: no emulator, the usual stopping rule", {
  rec <- seededRun(withPassing, mocks = emulatorMock)
  expect_null(rec$emu)
  expect_identical(rec$calls[[1]]$convergence, clusters::deoptimConvergence(NULL))
  expect_identical(unname(rec$calls[[1]]$DEoptimControl$initialpop[1:3, ]), unname(rec$passing))
})

test_that("without passing sets the emulator still starts the fit", {
  rec <- seededRun(list(thresh = 777, passingPars = NULL), params = list(initialPopMethod = "emulator"),
                   mocks = emulatorMock)
  expect_false(is.null(rec$emu))
  expect_identical(unname(rec$calls[[1]]$DEoptimControl$initialpop[1, ]), unname(rec$emuPop[1, ]))
})

test_that("more parameters than emulatorMaxParams, or a given initialpop, start the fit as 'calibration' does", {
  rec <- seededRun(withPassing, params = list(initialPopMethod = "emulator", emulatorMaxParams = 2L), mocks = emulatorMock)
  expect_null(rec$emu)
  expect_identical(rec$calls[[1]]$convergence, clusters::deoptimConvergence(NULL))
  expect_identical(unname(rec$calls[[1]]$DEoptimControl$initialpop[1:3, ]), unname(rec$passing))
  expect_message(expect_false(useEmulator(list(initialPopMethod = "emulator", emulatorMaxParams = 2L), 3L)),
                 "more than emulatorMaxParams")
  expect_true(useEmulator(list(initialPopMethod = "emulator", emulatorMaxParams = 18L), 3L))
  expect_false(useEmulator(list(initialPopMethod = "calibration", emulatorMaxParams = 18L), 3L))
  expect_error(useEmulator(list(initialPopMethod = "GP", emulatorMaxParams = 18L), 3L), "initialPopMethod")

  first <- seededRun(withPassing)
  given <- passingSets(first$lower, first$upper)[3:1, , drop = FALSE]
  rec <- seededRun(withPassing, params = list(initialPopMethod = "emulator", initialpop = given), mocks = emulatorMock)
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

test_that("the trust region is the box around the best fraction of the points, widened, scaled, inside the cube", {
  X <- rbind(c(0.40, 0.10), c(0.50, 0.30), c(0.90, 0.90), c(0.10, 0.95), c(0.45, 0.20),
             c(0.20, 0.80), c(0.70, 0.60), c(0.30, 0.50), c(0.60, 0.40), c(0.80, 0.70))
  score <- c(1, 2, 9, 8, 3, 7, 6, 5, 4, 10)
  box <- emulatorTrustBox(X, score, topFrac = 0.3, widen = 0.5)           # the 3 best: rows 1, 2, 5
  expect_equal(box, rbind(lower = c(0.35, 0), upper = c(0.55, 0.4)))     # x 0.4-0.5 +-0.05; y 0.1-0.3 +-0.1
  half <- emulatorTrustBox(X, score, topFrac = 0.3, widen = 0.5, scale = 0.5)
  expect_equal(half, rbind(lower = c(0.4, 0.1), upper = c(0.5, 0.3)))   # half the width about the same centre
  expect_equal(emulatorTrustBox(X, score, topFrac = 0.3, widen = 0)[, 1], c(lower = 0.4, upper = 0.5))
})

test_that("a trust round that improves by more than the noise shrinks the box; repeated failures stop the search", {
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
  ## the box never grows past the full trust box
  st <- emulatorTrustStep(list(scale = 0.9, successes = 0L, failures = 0L), improvement = 0, tol = 9, ctl)
  expect_equal(st$scale, 1)
  ## two successes needed before it shrinks
  ctl2 <- emulatorControlOf(list(trustSuccessTol = 2L))
  st <- emulatorTrustStep(list(scale = 1, successes = 0L, failures = 0L), improvement = 30, tol = 9, ctl2)
  expect_equal(st$scale, 1)
  expect_equal(emulatorTrustStep(st, improvement = 30, tol = 9, ctl2)$scale, ctl2$trustShrink)
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

test_that("the noise SD is pooled from the replicated points among the best", {
  U <- rbind(c(0.1, 0.1), c(0.1, 0.1), c(0.1, 0.1), c(0.2, 0.2), c(0.2, 0.2), c(0.9, 0.9), c(0.9, 0.9), c(0.5, 0.5),
             c(0.6, 0.6), c(0.7, 0.7))
  y <- c(100, 104, 102, 110, 106, 5000, 9000, 200, 300, 400)
  ## the 3 best distinct points (of 6) are (0.1, 0.1) (sd 2, 2 df) and (0.2, 0.2) (sd 2.83, 1 df), and (0.5, 0.5), not replicated; (0.9, 0.9) is not
  expect_equal(emulatorNoiseSD(U, y, topFrac = 0.34), sqrt((2 * 4 + 1 * 8) / 3))
  expect_true(is.na(emulatorNoiseSD(U[8:10, ], y[8:10], topFrac = 0.5)))  # no replicates
})

test_that("emulatorControl names are checked", {
  expect_identical(emulatorControlOf(list(maxEvals = 10L))$maxEvals, 10L)
  expect_identical(emulatorControlOf()$designPerParam, emulatorDefaults()$designPerParam)
  expect_error(emulatorControlOf(list(maxEval = 10L)), "unknown element")
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

test_that("the local phase is trust rounds that shrink around the minimum and stop when they improve less than the noise", {
  ## one global batch: the rounds find the minimum, each improving more than the noise and halving the region
  res <- quadEmulate(utils::modifyList(quadCtl, list(globalMaxEvals = 76L, maxEvals = 2000L)))
  st <- res$status[res$status$phase == "local"]
  expect_gte(max(st$round), 3L)
  expect_lt(res$nEvaluated, 2000L)                                       # stopped by the rule, not the budget
  expect_lt(utils::tail(st$boxWidth, 1L), st$boxWidth[1L] / 4)
  expect_lt(min(st$scale), 0.5)
  expect_lt(min(quadValue(res$pop)) - 100, 5)
  ## two global batches already find it: the rounds improve less than the noise, and the search stops
  res <- quadEmulate(utils::modifyList(quadCtl, list(globalMaxEvals = 88L, maxEvals = 2000L)))
  expect_identical(res$stopReason, "improvement below noise")
  expect_identical(max(res$status$round, na.rm = TRUE), 2L)
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
