## The DEoptim cluster's size is also its population size: clusters:::.clusterNP() sets NP to the
## number of workers the cluster was built with, discarding any NP the caller asked for. So the only
## way to choose NP is to choose the cluster size, which is `nCoresNeeded` in
## fireSenseUtils::runDEoptim() -- and this module never passed it, so the package default
## (10 * length(lower), i.e. 120 for a 12-parameter fit) always won.
##
## FireSense phase 2 (2026-09-16): measured per-ELF, a generation costs the SLOWEST of NP evaluations,
## and the expected maximum barely moves with NP (68.0 s at NP 120 vs 64.9 s at NP 60 for ELF 4.1).
## Throughput therefore comes from running more ELFs at once, not from shorter generations: at 704
## cores, NP 60 allows 11 concurrent fits where NP 120 allows 5. Choosing NP needs it to be a
## parameter of this module.

test_that("the module declares an nCoresNeeded parameter", {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_spreadFit.R"), keep.source = FALSE)
  code <- paste(vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1)),
                collapse = "\n")
  expect_match(code, "defineParameter\\(\"nCoresNeeded\"")
})

test_that("the runDEoptim() call passes nCoresNeeded from the parameter", {
  for (nc in list(30L, NULL)) {
    sim <- toySim(list(nCoresNeeded = nc))
    rec <- mockFitAndLedger(sim)
    suppressMessages(SpaDES.core::spades(sim))
    expect_identical(rec$deArgs$nCoresNeeded, nc)
  }
})

## 2026-11-17 (GPSC, null arm, ELFs 14.3 and 11.3): `cores = 40L` with 3 parameters. clusterSetup() builds 40 workers
## and sets NP to 40, but the module's population size was deoptimPopulationSize(3) = 30, so the emulator and the start
## population had 30 rows: "initialpop has 30 rows but the population (the number of workers) is 40".
test_that("a single number of cores is the population size: the emulator and the start population have that many rows", {
  rec <- seededRun(params = list(cores = 40L, initialPopMethod = "emulator"), mocks = emulatorMock)
  expect_identical(rec$emu$NP, 40L)
  expect_identical(nrow(rec$emuPop), 40L)
  expect_gte(nrow(rec$calls[[1]]$DEoptimControl$initialpop), 40L)
})

test_that("the population size is nCoresNeeded (else the default) unless cores is a number of workers", {
  np <- function(...) spreadFitNP(toySim(list(...)))
  ## the module default, cores = 1L, is one host: clusterSetup() sizes it by nCoresNeeded
  lower <- c(a = 0, b = 0, c = 0, d = 0)
  expect_identical(np(lower = lower), fireSenseUtils::deoptimPopulationSize(4L))
  expect_identical(np(cores = 1L, nCoresNeeded = 25L), 25L)
  ## host names: nCoresNeeded, as before
  expect_identical(np(cores = c("hostA", "hostB"), nCoresNeeded = 25L), 25L)
  ## a number of workers is the cluster's size, whatever nCoresNeeded says
  expect_identical(np(cores = 40L, nCoresNeeded = 25L), 40L)
  expect_identical(np(cores = 12L), 12L)
  ## a host named more than once is that many workers (clusters::requestedWorkers())
  expect_identical(np(cores = rep("hostA", 40), nCoresNeeded = 25L), 40L)
  expect_identical(np(cores = c(rep("hostA", 3), "hostB")), 4L)
  ## no cluster: nCoresNeeded or the default
  expect_identical(np(cores = NULL, nCoresNeeded = 25L), 25L)
  expect_identical(np(cores = NA, nCoresNeeded = 25L), 25L)
})

## 2026-11-18: cores = rep("hostA", 40) with 3 parameters: clusterSetup() asks for 40 workers (a repeated host is that
## many), but spreadFitNP() gave deoptimPopulationSize(3) = 30 rows
test_that("a repeated host name is the population size: the emulator and the start population have that many rows", {
  rec <- seededRun(params = list(cores = rep("hostA", 40), initialPopMethod = "emulator"), mocks = emulatorMock)
  expect_identical(rec$emu$NP, 40L)
  expect_identical(nrow(rec$emuPop), 40L)
  expect_gte(nrow(rec$calls[[1]]$DEoptimControl$initialpop), 40L)
})
