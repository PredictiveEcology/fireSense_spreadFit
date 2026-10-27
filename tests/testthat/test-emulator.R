## initialPopMethod = "emulator": a Gaussian-process emulator of the objective chooses DEoptim's starting population
## (prototype 2026-10-25/26: ELF 6.3.1 seeded reached in 16 generations the median the unseeded fit took 249 for).

## The emulator's search, replaced by a recorder that returns NP points at 0.3 of the way from lower to upper
emulatorMock <- function(rec) list(
  emulatorStartPopulation = function(objArgs, lower, upper, passingPars, NP, control, seed, cores, nCoresNeeded, ...) {
    rec$emu <- list(objArgs = objArgs, passingPars = passingPars, NP = NP, cores = cores, nCoresNeeded = nCoresNeeded)
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
  expect_identical(call$convergence, clusters::deoptimConvergence(list(window = 100L, minGenerations = 100L)))
  ## the box is centred on the calibration's passing sets; the population size and the cores are the fit's
  expect_identical(rec$emu$passingPars, rec$passing)
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
  expect_identical(rec$calls[[1]]$convergence, clusters::deoptimConvergence(list(window = 250L, minGenerations = 100L)))
})

test_that("initialPopMethod = 'calibration' (the default) is the fit as before: no emulator, the usual stopping rule", {
  rec <- seededRun(withPassing, mocks = emulatorMock)
  expect_null(rec$emu)
  expect_identical(rec$calls[[1]]$convergence, clusters::deoptimConvergence(NULL))
  expect_identical(unname(rec$calls[[1]]$DEoptimControl$initialpop[1:3, ]), unname(rec$passing))
})

test_that("without passing sets the emulator searches between lower and upper", {
  rec <- seededRun(list(thresh = 777, passingPars = NULL), params = list(initialPopMethod = "emulator"),
                   mocks = emulatorMock)
  expect_null(rec$emu$passingPars)
  expect_identical(emulatorBox(NULL, rec$lower, rec$upper), rbind(lower = rec$lower, upper = rec$upper))
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

test_that("emulatorBox() is the passing sets' range widened, at least a minimum width, inside the bounds", {
  lo <- c(a = 0, b = -10, c = 0); up <- c(a = 1, b = 10, c = 100)
  passing <- rbind(c(a = 0.4, b = 5, c = 50), c(a = 0.5, b = 9, c = 50))
  box <- emulatorBox(passing, lo, up, widen = 1, minWidthFrac = 0.1)
  expect_identical(colnames(box), names(lo))
  expect_equal(box[, "a"], c(lower = 0.3, upper = 0.6))     # range 0.1, widened by 0.1 each side
  expect_equal(box[, "b"], c(lower = 1, upper = 10))        # 3 to 11, clipped to the bound 10
  expect_equal(box[, "c"], c(lower = 45, upper = 55))       # one value: the minimum width, 10% of 100
})

test_that("emulatorControl names are checked", {
  expect_identical(emulatorControlOf(list(maxEvals = 10L))$maxEvals, 10L)
  expect_identical(emulatorControlOf()$designPerParam, emulatorDefaults()$designPerParam)
  expect_error(emulatorControlOf(list(maxEval = 10L)), "unknown element")
})

test_that("the emulator finds points near the minimum of a noisy quadratic within its budget", {
  lower <- c(a = 0, b = -5, c = 0, d = -10); upper <- c(a = 1, b = 5, c = 20, d = 0)
  opt <- c(a = 0.3, b = 1, c = 12, d = -2)
  trueValue <- function(P)
    100 + 1000 * rowSums(sweep(sweep(P[, names(lower), drop = FALSE], 2, opt), 2, upper - lower, "/")^2)
  nEvaluated <- 0L
  noisy <- function(P, seed) {
    expect_identical(colnames(P), names(lower))
    nEvaluated <<- nEvaluated + nrow(P)
    withr::with_seed(seed, trueValue(P) + stats::rnorm(nrow(P), 0, 2))
  }
  ctl <- list(designPerParam = 10L, globalMaxEvals = 150L, maxEvals = 250L, minEvals = 100L, batchSize = 12L,
              nCandidates = 500L, nPerturb = 300L, gpThreads = 1L)
  res <- suppressMessages(emulateObjective(noisy, lower, upper, NP = 8L, control = ctl, seed = 1L))
  expect_identical(colnames(res$pop), names(lower))
  expect_identical(nrow(res$pop), 8L)
  expect_false(anyDuplicated(res$pop) > 0)
  expect_true(all(t(res$pop) >= lower & t(res$pop) <= upper))
  expect_lte(nEvaluated, ctl$maxEvals + ctl$batchSize)
  expect_identical(nrow(res$evals), nEvaluated)
  expect_true(all(c("global", "local") %in% res$status$phase))
  ## the start of the box is 100-1500 above the minimum; the best points are within a few units of it
  expect_lt(median(trueValue(res$pop)) - 100, 10)
  expect_identical(res$objFunVal, sort(res$objFunVal))
})

test_that("the emulator's evaluations run on a cluster built by clusters::clusterSetup() from the fit's cores", {
  rec <- new.env()
  testthat::local_mocked_bindings(.package = "clusters", clusterSetup = function(...) {
    a <- list(...)
    rec$setup <- a
    rec$shipped <- mget(unlist(a$objsNeeded), envir = a$envir)
    list(cluster = list("w1", "w2", "w3"))
  })
  testthat::local_mocked_bindings(emulateObjective = function(evaluate, lower, upper, box, NP, control, seed) {
    rec$batchSize <- control$batchSize
    list(pop = NULL)
  })
  objArgs <- list(landscape = "L", annualDTx1000 = 1, nonAnnualDTx1000 = 2, fireBufferedListDT = 3,
                  historicalFires = 4, mutuallyExclusive = list(), tests = "adTest")
  lo <- c(a = 0); up <- c(a = 1)
  suppressMessages(emulatorStartPopulation(objArgs, lo, up, NULL, NP = 30L, control = list(), seed = 1L,
                                           cores = c("hostA", "hostB"), nCoresNeeded = 30L, runName = "toy"))
  expect_identical(rec$setup$cores, c("hostA", "hostB"))
  expect_identical(rec$setup$nCoresNeeded, 30L)
  expect_identical(rec$shipped, objArgs[emulatorDataArgs])
  expect_identical(rec$batchSize, 3L)                                   # one evaluation per worker per batch
  ## without cores, nothing is built and the evaluations run here
  rec$setup <- NULL
  suppressMessages(emulatorStartPopulation(objArgs, lo, up, NULL, NP = 30L, control = list(), seed = 1L, cores = NA))
  expect_null(rec$setup)
  expect_null(rec$batchSize)
})
