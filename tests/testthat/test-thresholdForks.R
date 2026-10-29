source(testthat::test_path("..", "..", "R", "thresholdForks.R"))

test_that("a large heap on a shared host yields as many forks as memory allows, no floor of four", {
  # 2026-09-07: ~90 GB parent, ~300 GB left on a 1 TB host shared with two other
  # simulations. The old code returned 4 here (its floor); 4 x 90 GB x 3 sims
  # exhausted the host.
  expect_equal(thresholdForks(heapMB = 90000, availMB = 300000, nPars = 96,
                              detCores = 68, activeThreads = 0L), 3L)
})

test_that("options(mc.cores) is an explicit cap and is honoured", {
  expect_equal(thresholdForks(heapMB = 90000, availMB = 900000, nPars = 96,
                              detCores = 68, activeThreads = 0L, mcCores = 2L), 2L)
})

test_that("busy threads reduce the CPU budget only, not the memory budget", {
  # 68 cores, 20 busy: 34 - 20 = 14 forks; the old code subtracted 20 from
  # min(everything) and fell back on its floor.
  expect_equal(thresholdForks(heapMB = 1000, availMB = 1e6, nPars = 96,
                              detCores = 68, activeThreads = 20L), 14L)
})

test_that("a small heap on an idle host uses half the cores", {
  expect_equal(thresholdForks(heapMB = 20, availMB = 1e6, nPars = 96, detCores = 68), 34L)
})

test_that("the number of parameter sets caps the forks", {
  expect_equal(thresholdForks(heapMB = 20, availMB = 1e6, nPars = 2, detCores = 68), 2L)
})

test_that("never fewer than one fork; unknown memory means one", {
  expect_equal(thresholdForks(heapMB = 90000, availMB = 300000, nPars = 96,
                              detCores = 68, activeThreads = 100L), 1L)
  expect_equal(thresholdForks(heapMB = 20, availMB = NULL, nPars = 96, detCores = 68), 1L)
  expect_equal(thresholdForks(heapMB = 20, availMB = numeric(0), nPars = 96, detCores = 68), 1L)
  expect_equal(thresholdForks(heapMB = 20, availMB = NA_real_, nPars = 96, detCores = 68), 1L)
})

test_that("forkBudget() honours options(mc.cores), whoever calls it", {
  withr::local_options(mc.cores = 1L)
  expect_equal(forkBudget(10000L), 1L)
  withr::local_options(mc.cores = 3L)
  expect_lte(forkBudget(10000L), 3L)
  expect_equal(forkBudget(1L), 1L)
})

test_that("loadCap lowers the fork count; the default Inf leaves it unchanged", {
  expect_equal(thresholdForks(heapMB = 20, availMB = 1e6, nPars = 96, detCores = 68, loadCap = 5L), 5L)
  expect_equal(thresholdForks(heapMB = 20, availMB = 1e6, nPars = 96, detCores = 68, loadCap = Inf), 34L)
  expect_equal(thresholdForks(heapMB = 20, availMB = 1e6, nPars = 96, detCores = 68, loadCap = 100L), 34L)
})

test_that("localSingleThreadBLAS() gives one BLAS thread until the caller returns", {
  skip_if_not_installed("RhpcBLASctl")
  before <- RhpcBLASctl::blas_get_num_procs()
  inside <- (function() {
    localSingleThreadBLAS()
    RhpcBLASctl::blas_get_num_procs()
  })()
  expect_equal(inside, 1L)
  expect_equal(RhpcBLASctl::blas_get_num_procs(), before)
})

test_that("forkBudget() sizes the forks by the CPUs the job was given, not the node's", {
  # 2026-10-28: a Slurm job of part of a 64-CPU node forked 32 children (half the node's cores).
  # parallelly reads the environment once per session, so ask a fresh R process that starts inside the "job".
  skip_if_not_installed("callr")
  skip_if(parallel::detectCores() <= 2L, "needs a machine with more than 2 cores")
  forks <- callr::r(function(src) {
    options(mc.cores = NULL)
    source(src)
    forkBudget(10000L)
  }, args = list(testthat::test_path("..", "..", "R", "thresholdForks.R")),
  env = c(callr::rcmd_safe_env(), SLURM_JOB_ID = "1", SLURM_CPUS_ON_NODE = "2", SLURM_CPUS_PER_TASK = "2"))
  expect_lte(forks, 2L)
})

test_that("forkBudget() does not exceed the memory the job was given", {
  withr::local_options(mc.cores = NULL)
  heapMB <- sum(gc()[, 2])
  testthat::local_mocked_bindings(.slurmAllocation = function(...) list(memLimitGB = 1, memAvailableGB = 1),
                                  .package = "clusters")
  expect_lte(forkBudget(10000L), max(1L, floor(1024 / heapMB)))
})

test_that("forkBudget() ignores the node's busy threads and load in a Slurm job given part of the node", {
  # 2026-10-29, GPSC: 42 CPUs of a 64-CPU node shared with other jobs logged "Using 3 cores" for the threshold
  # calibration. Other jobs' threads run on CPUs the job's cgroup keeps away from it.
  withr::local_options(mc.cores = NULL)
  testthat::local_mocked_bindings(.slurmAllocation = function(...) list(memLimitGB = 1e6, memAvailableGB = 1e6),
                                  numActiveThreads = function(...) 100L, .package = "clusters")
  testthat::local_mocked_bindings(nThreads = function(...) 1L, .package = "fireSenseUtils")
  testthat::local_mocked_bindings(availableCores = function(...) 8L, .package = "parallelly")
  testthat::local_mocked_bindings(detectCores = function(...) 64L, .package = "parallel")
  expect_equal(forkBudget(10000L), 4L)  # half of the job's 8 CPUs
  # a whole-node job still sees the node's load
  testthat::local_mocked_bindings(availableCores = function(...) 64L, .package = "parallelly")
  expect_equal(forkBudget(10000L), 1L)
})

test_that("forkBudget() outside Slurm still subtracts busy threads and caps by load", {
  withr::local_options(mc.cores = NULL)
  testthat::local_mocked_bindings(.slurmAllocation = function(...) NULL, numActiveThreads = function(...) 2L,
                                  .package = "clusters")
  testthat::local_mocked_bindings(nThreads = function(...) 100L, .package = "fireSenseUtils")
  testthat::local_mocked_bindings(availableCores = function(...) 8L, .package = "parallelly")
  testthat::local_mocked_bindings(detectCores = function(...) 64L, .package = "parallel")
  expect_equal(forkBudget(10000L), 2L)  # floor(8 * 0.5) - 2 busy
  testthat::local_mocked_bindings(numActiveThreads = function(...) 0L, .package = "clusters")
  testthat::local_mocked_bindings(nThreads = function(...) 1L, .package = "fireSenseUtils")
  expect_equal(forkBudget(10000L), 1L)  # load cap
})
