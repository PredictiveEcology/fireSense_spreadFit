#' How many forks the threshold self-calibration may use
#'
#' `mcmapply()` forks this process. Each child starts as a copy-on-write image
#' of the parent and diverges as R's garbage collector marks the heap, so the
#' honest per-child budget is the parent's own R heap, not a fixed 5000 MB: on
#' 2026-09-07 the forks of a ~90 GB parent each reached ~37 GB of private
#' memory within minutes. There is no floor above one. The previous
#' `max(4L, ...)` forced four such children per simulation whatever the host
#' had left and whatever `options(mc.cores)` said; three simulations doing that
#' on one host exhausted 1 TB.
#'
#' Threads already busy on the host reduce the CPU budget only. The previous
#' code subtracted them from the minimum of every budget, which on a shared
#' host went negative and landed back on the floor.
#'
#' @param heapMB numeric; R heap in use in this process, MB (`sum(gc()[, 2])`).
#' @param availMB numeric or `NULL`; memory available on the host, MB.
#'   Unknown (`NULL`, empty, non-finite) means one fork.
#' @param nPars integer; number of parameter sets to evaluate.
#' @param detCores integer; `parallel::detectCores()`.
#' @param activeThreads integer; busy threads already running on the host.
#' @param mcCores integer or `NULL`; `getOption("mc.cores")`, an explicit cap.
#' @param loadCap numeric; cap from the host's load (`fireSenseUtils::nThreads()`). `Inf` (default) is no cap.
#' @return integer, at least 1.
thresholdForks <- function(heapMB, availMB, nPars, detCores, activeThreads = 0L,
                           mcCores = NULL, loadCap = Inf) {
  byMem <- if (!length(availMB) || !is.finite(availMB[1])) 1L
           else floor(availMB[1] / max(heapMB, 1))
  byCPU <- floor(detCores * 0.5) - activeThreads
  n <- min(c(byMem, byCPU, nPars, mcCores, loadCap))
  max(1L, as.integer(n))
}

#' How many forks this process may start now
#'
#' [thresholdForks()] with this host's state: the R heap in use (each fork is a copy-on-write image of this
#' process that diverges as the garbage collector marks the heap, so budget one heap per fork; the `gc()` call
#' also leaves less garbage to copy), the memory available, the busy threads, `detectCores()` and
#' `options(mc.cores)`, and capped by the host's load (`fireSenseUtils::nThreads()`). Used by the threshold self-calibration and the `"powerLaw"` calibration, so several
#' fits on one host, or a fit run inside a cluster worker, share the same caps.
#'
#' @param nTasks integer; number of tasks to spread over the forks.
#' @return integer, at least 1.
forkBudget <- function(nTasks) {
  heapMB <- sum(gc()[, 2])
  availMB <- pemisc::availableMemory()
  availMB <- if (length(availMB)) availMB / 1e6 else NULL
  thresholdForks(heapMB = heapMB, availMB = availMB, nPars = nTasks, detCores = parallel::detectCores(),
                 activeThreads = clusters::numActiveThreads(), mcCores = getOption("mc.cores"),
                 loadCap = fireSenseUtils::nThreads())
}

#' One BLAS thread until the calling function returns
#'
#' A fork inherits the parent's OpenBLAS thread pool, so each fork of the threshold calibration ran several
#' threads (64 on an 80-core host) and five fits calibrating at once pushed the 1-minute load past 370.
#' Call this before forking: the forks inherit one thread, and the previous number is restored when the
#' caller exits.
#'
#' @param envir environment whose exit restores the previous number of BLAS threads.
#' @param threads integer; BLAS threads until then (the emulator's GP fits take a few, `emulatorDefaults()`).
#' @return `NULL`, invisibly.
localSingleThreadBLAS <- function(envir = parent.frame(), threads = 1L) {
  previous <- RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(threads)
  withr::defer(RhpcBLASctl::blas_set_num_threads(previous), envir = envir)
  invisible(NULL)
}
