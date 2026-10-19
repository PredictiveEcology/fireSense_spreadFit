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
#' @return integer, at least 1.
thresholdForks <- function(heapMB, availMB, nPars, detCores, activeThreads = 0L,
                           mcCores = NULL) {
  byMem <- if (!length(availMB) || !is.finite(availMB[1])) 1L
           else floor(availMB[1] / max(heapMB, 1))
  byCPU <- floor(detCores * 0.5) - activeThreads
  n <- min(c(byMem, byCPU, nPars, mcCores))
  max(1L, as.integer(n))
}

#' How many forks this process may start now
#'
#' [thresholdForks()] with this host's state: the R heap in use (each fork is a copy-on-write image of this
#' process that diverges as the garbage collector marks the heap, so budget one heap per fork; the `gc()` call
#' also leaves less garbage to copy), the memory available, the busy threads, `detectCores()` and
#' `options(mc.cores)`. Used by the threshold self-calibration and the `"powerLaw"` calibration, so several
#' fits on one host, or a fit run inside a cluster worker, share the same caps.
#'
#' @param nTasks integer; number of tasks to spread over the forks.
#' @return integer, at least 1.
forkBudget <- function(nTasks) {
  heapMB <- sum(gc()[, 2])
  availMB <- pemisc::availableMemory()
  availMB <- if (length(availMB)) availMB / 1e6 else NULL
  thresholdForks(heapMB = heapMB, availMB = availMB, nPars = nTasks, detCores = parallel::detectCores(),
                 activeThreads = clusters::numActiveThreads(), mcCores = getOption("mc.cores"))
}
