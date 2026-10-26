## A Gaussian-process emulator of the spread objective, used to choose DEoptim's starting population
## (`initialPopMethod = "emulator"`). Ported from the 2026-10-25/26 prototype (ELFs 6.6.1, 6.6.2, 6.3.1, 6.1.2):
## ~1,000 real evaluations found sets whose DEoptim population then reached, within 16 generations, a median the
## unseeded fit needed 249 generations for (6.3.1). The emulator knows nothing about the covariates or the link: it
## sees only `lower`/`upper` (named) and a function that evaluates parameter sets.

#' The emulator's settings and their defaults
#'
#' `emulatorControl` (a module parameter) overrides any of these by name.
#'
#' @return a named list:
#'   - `designPerParam` (20): the space-filling design (maximin Latin hypercube, `lhs::maximinLHS()`) has this many
#'     points per parameter, each evaluated once;
#'   - `designRepPoints` (4), `designReps` (6): that many of the design points are evaluated `designReps` times
#'     more, so the GP sees the noise away from the optimum;
#'   - `globalMaxEvals` (1000), `maxEvals` (1600): the budget of real evaluations of the global phase, and of both
#'     phases together;
#'   - `minEvals` (600): no phase stops for lack of improvement before this many evaluations;
#'   - `globalWindow` (6), `localWindow` (5), `relTol` (0.002): a phase stops when the GP's best mean (at an
#'     evaluated point) improved by less than `relTol` times itself over its last `globalWindow`/`localWindow`
#'     batches (the prototype: 2 objective units at ~938);
#'   - `batchSize` (`NULL`): evaluations per batch; `NULL` is the cluster's number of workers, or 48 without one;
#'   - `floorFrac` (0.8), `globalCapRatio` (10), `localCapRatio` (1.18): the GP models
#'     `log(min(y, cap) - floorFrac * min(y))`, with `cap` `globalCapRatio` times the best value seen in the global
#'     phase and `localCapRatio` times the GP's best mean in the local phase, so failures (the objective's fail
#'     value, and runaway sets scored in the millions) sit on a plateau instead of dominating the fit;
#'   - `nCandidates` (4000), `nPerturb` (3000), `perturbSD` (c(0.005, 0.015, 0.04, 0.1)), `nTop` (20): each batch
#'     is chosen from `nCandidates` Latin-hypercube points (in the trust region in the local phase) and `nPerturb`
#'     Gaussian perturbations, at each `perturbSD` (unit-box scale), of the `nTop` best points by GP mean;
#'   - `exclusion` (0.25): batch members are at least this far apart in lengthscale units;
#'   - `lcbSD` (2): the exploration picks minimise the GP mean minus this many standard deviations;
#'   - `trustTop` (40): the local phase's trust region is the box around the `trustTop` best points by GP mean,
#'     widened by half its width each side;
#'   - `boxWiden` (1), `boxMinWidthFrac` (0.1): the search box is the range of the threshold calibration's passing
#'     sets widened by `boxWiden` times that range each side, at least `boxMinWidthFrac` of `upper - lower` wide,
#'     clipped to `lower`/`upper`;
#'   - `thetaBounds` (c(0.01, 5)), `gpMaxit` (200), `gpThreads` (6): the GP's lengthscale bounds (unit box),
#'     optimiser iterations (`hetGP::mleHetGP()`) and BLAS threads (the rest of the search runs on one).
emulatorDefaults <- function() {
  list(designPerParam = 20L, designRepPoints = 4L, designReps = 6L,
       globalMaxEvals = 1000L, maxEvals = 1600L, minEvals = 600L,
       globalWindow = 6L, localWindow = 5L, relTol = 0.002, batchSize = NULL,
       floorFrac = 0.8, globalCapRatio = 10, localCapRatio = 1.18,
       nCandidates = 4000L, nPerturb = 3000L, perturbSD = c(0.005, 0.015, 0.04, 0.1), nTop = 20L,
       exclusion = 0.25, lcbSD = 2, trustTop = 40L, boxWiden = 1, boxMinWidthFrac = 0.1,
       thetaBounds = c(0.01, 5), gpMaxit = 200L, gpThreads = 6L)
}

## `control` merged onto the defaults; an unknown name is an error, not a silently ignored setting
emulatorControlOf <- function(control = list()) {
  ctl <- emulatorDefaults()
  unknown <- setdiff(names(control), names(ctl))
  if (length(unknown))
    stop("emulatorControl has unknown element(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  ctl[names(control)] <- control
  ctl
}

#' The box the emulator searches
#'
#' Centred on the threshold calibration's passing sets: their range, widened by `widen` times itself each side and
#' at least `minWidthFrac` of the bounds' width, clipped to `lower`/`upper`. Without passing sets, `lower`/`upper`.
#'
#' @param passingPars matrix with the columns of `lower`, or `NULL`.
#' @param lower,upper the bounds, named.
#' @param widen,minWidthFrac see `emulatorDefaults()` (`boxWiden`, `boxMinWidthFrac`).
#' @return a 2-row matrix (`lower`, `upper`), columns `names(lower)`.
emulatorBox <- function(passingPars, lower, upper, widen = 1, minWidthFrac = 0.1) {
  box <- rbind(lower = lower, upper = upper)
  if (is.null(passingPars) || !NROW(passingPars)) return(box)
  r <- apply(as.matrix(passingPars)[, names(lower), drop = FALSE], 2L, range)
  w <- pmax((r[2, ] - r[1, ]) * (1 + 2 * widen), minWidthFrac * (upper - lower))
  mid <- (r[1, ] + r[2, ]) / 2
  box[1, ] <- pmax(lower, mid - w / 2)
  box[2, ] <- pmin(upper, mid + w / 2)
  box
}

#' Search a noisy objective with a Gaussian-process emulator and return its best evaluated points
#'
#' A space-filling design (`lhs::maximinLHS()`, plus replicates) in `box`, then batches chosen on a GP
#' (`hetGP::mleHetGP()`, Matern 5/2) fitted to every evaluation so far, on a capped log scale (see
#' `emulatorDefaults()`): expected improvement over the GP's best mean (greedy, with an exclusion distance), a lower
#' confidence bound for exploration, and replicates of the best points. A global phase over the box, then a local
#' phase in a trust region around the best points, each until the GP's best mean stops improving or the budget is
#' spent.
#'
#' @param evaluate a function of a matrix of parameter sets (one per row, columns `names(lower)`) and an integer
#'   `seed`, returning one objective value per row (`emulatorEvaluator()`).
#' @param lower,upper the bounds, named.
#' @param box the search box, `emulatorBox()`.
#' @param NP integer; how many points to return.
#' @param control list; see `emulatorDefaults()`.
#' @param seed integer; the design, the candidates and the evaluations' seeds follow from it.
#' @return a list: `pop`, the `NP` best distinct evaluated points by their mean value over their evaluations (fewer
#'   if fewer were evaluated), columns `names(lower)`; `objFunVal`, those means; `evals`, every evaluation (`batch`,
#'   `kind`, `value` and the parameters); `status`, one row per batch.
emulateObjective <- function(evaluate, lower, upper, box = rbind(lower = lower, upper = upper), NP,
                             control = list(), seed = 1L) {
  ctl <- emulatorControlOf(control)
  parNames <- names(lower)
  d <- length(lower)
  lo <- box[1L, parNames]; width <- box[2L, parNames] - lo
  fromUnit <- function(U) {
    P <- sweep(sweep(as.matrix(U), 2L, width, "*"), 2L, lo, "+")
    colnames(P) <- parNames
    P
  }
  B <- if (is.null(ctl$batchSize)) 48L else as.integer(ctl$batchSize)
  localSingleThreadBLAS() # the evaluations, and the GP's predictions; its fits take `gpThreads`
  withr::with_seed(seed, {
    U <- lhs::maximinLHS(as.integer(ctl$designPerParam * d), d)
    repRows <- unique(round(seq(1, nrow(U), length.out = ctl$designRepPoints + 1L)))[seq_len(ctl$designRepPoints)]
    U <- rbind(U, U[rep(repRows, each = ctl$designReps), , drop = FALSE])
    kind <- c(rep("design", nrow(U) - length(repRows) * ctl$designReps), rep("rep", length(repRows) * ctl$designReps))
    y <- evaluate(fromUnit(U), seed = seed)
    batch <- rep(0L, nrow(U))
    status <- NULL
    phase <- "global"
    k <- 0L
    repeat {
      if (nrow(U) >= ctl$maxEvals) break
      if (phase == "global" && (nrow(U) >= ctl$globalMaxEvals || emulatorConverged(status, phase, nrow(U), ctl)))
        phase <- "local"
      if (phase == "local" && emulatorConverged(status, phase, nrow(U), ctl)) break
      k <- k + 1L
      local <- phase == "local"
      best <- if (is.null(status)) min(y, na.rm = TRUE) else min(status$bestMean)
      cap <- if (local) ctl$localCapRatio * best else ctl$globalCapRatio * min(y, na.rm = TRUE)
      sc <- emulatorScale(y, cap, ctl$floorFrac)
      m <- fitEmulatorGP(U, sc$z, ctl)
      sel <- selectEmulatorBatch(m, U, B, local, ctl)
      new <- evaluate(fromUnit(sel$U), seed = seed + 1e6 * k)
      U <- rbind(U, sel$U); y <- c(y, new); kind <- c(kind, sel$kind); batch <- c(batch, rep(k, nrow(sel$U)))
      status <- rbind(status, data.table(batch = k, phase = phase, nEval = nrow(U),
                                         bestMean = sc$untr(sel$bestZ), batchBest = min(new, na.rm = TRUE)))
      message("emulator batch ", k, " (", phase, "): ", nrow(U), " evaluations; GP best mean ",
              signif(utils::tail(status$bestMean, 1L), 6), "; best single value ", signif(min(y, na.rm = TRUE), 6))
    }
  })
  P <- fromUnit(U)
  scores <- data.table(member = seq_along(y), rep = 1L, value = y)[!is.na(value)]
  best <- fireSenseUtils::bestByReplicatedMean(P, scores, n = NP)
  pop <- as.matrix(best$params)
  colnames(pop) <- parNames
  list(pop = pop, objFunVal = best$objFunVal,
       evals = data.table(batch = batch, kind = kind, value = y, P), status = status)
}

## The GP's response: log(min(y, cap) - floor), floor `floorFrac` of the best value. A missing value (an evaluation
## that returned none) is at the cap. `untr()` maps a GP mean back to objective units.
emulatorScale <- function(y, cap, floorFrac = 0.8) {
  best <- min(y, na.rm = TRUE)
  if (!is.finite(best) || best <= 0)
    stop("the emulator models the log of a positive objective; the best value is ", best, call. = FALSE)
  floor <- floorFrac * best
  cap <- max(cap, best)
  y[is.na(y)] <- cap
  list(z = log(pmin(y, cap) - floor), untr = function(z) floor + exp(z))
}

## The GP fitted to every evaluation, on the unit box (hetGP pools replicates)
fitEmulatorGP <- function(U, z, ctl) {
  localSingleThreadBLAS(threads = ctl$gpThreads)
  d <- ncol(U)
  hetGP::mleHetGP(X = U, Z = z, lower = rep(ctl$thetaBounds[1L], d), upper = rep(ctl$thetaBounds[2L], d),
                  covtype = "Matern5_2", maxit = ctl$gpMaxit, settings = list(trace = -1))
}

## Whether `phase` has stopped improving: at least `minEvals` evaluations, more than `window` batches in the phase,
## and the GP's best mean improved by less than `relTol` times itself over the last `window` of them
emulatorConverged <- function(status, phase, nEval, ctl) {
  if (is.null(status) || nEval < ctl$minEvals) return(FALSE)
  b <- status$bestMean[status$phase == phase]
  window <- if (phase == "global") ctl$globalWindow else ctl$localWindow
  n <- length(b)
  n > window && b[n - window] - b[n] < ctl$relTol * b[n]
}

#' One batch of points to evaluate, chosen on the GP
#'
#' Expected improvement over the GP's best mean at an evaluated point (greedy, each pick at least
#' `exclusion` lengthscales from the others), then a lower confidence bound (exploration), then replicates of the
#' best evaluated points (fewest replicates first). In the local phase the Latin-hypercube candidates are drawn in
#' the trust region (see `emulatorDefaults()`) and more replicates are taken.
#'
#' @param m the GP (`fitEmulatorGP()`).
#' @param U every evaluated point so far, unit box.
#' @param B integer; batch size.
#' @param local logical; the local phase.
#' @param ctl `emulatorControlOf()`.
#' @return list: `U` (the batch, unit box), `kind` ("ei", "lcb", "rep"), `bestZ` (the GP's best mean at an
#'   evaluated point, response scale).
selectEmulatorBatch <- function(m, U, B, local, ctl) {
  nRep <- round(B / if (local) 4 else 6); nLCB <- round(B / 8); nEI <- B - nRep - nLCB
  d <- ncol(U)
  keys <- apply(U, 1L, function(r) paste(sprintf("%a", r), collapse = ","))
  X0 <- U[!duplicated(keys), , drop = FALSE]
  mult <- tabulate(match(keys, keys[!duplicated(keys)]), nrow(X0))
  pm <- predict(m, X0)$mean
  ord <- order(pm)
  top <- X0[ord[seq_len(min(ctl$nTop, nrow(X0)))], , drop = FALSE]
  cand <- lhs::randomLHS(ctl$nCandidates, d)
  if (local) {
    tb <- X0[ord[seq_len(min(ctl$trustTop, nrow(X0)))], , drop = FALSE]
    tlo <- apply(tb, 2L, min); thi <- apply(tb, 2L, max); w <- thi - tlo
    tlo <- pmax(tlo - w / 2, 0); thi <- pmin(thi + w / 2, 1)
    cand <- sweep(sweep(cand, 2L, thi - tlo, "*"), 2L, tlo, "+")
  }
  cand <- rbind(cand, do.call(rbind, lapply(ctl$perturbSD, function(s) {
    base <- top[sample(nrow(top), ctl$nPerturb, replace = TRUE), , drop = FALSE]
    pmin(pmax(base + matrix(stats::rnorm(length(base), 0, s), ncol = d), 0), 1)
  })))
  p <- predict(m, cand)
  s <- sqrt(pmax(p$sd2, 1e-12))
  fstar <- min(pm)
  u <- (fstar - p$mean) / s
  ei <- (fstar - p$mean) * stats::pnorm(u) + s * stats::dnorm(u)
  lcb <- p$mean - ctl$lcbSD * s
  theta <- rep_len(m$theta, d)
  far <- function(A, b) sqrt(colSums(((t(A) - b) / theta)^2)) > ctl$exclusion
  greedy <- function(score, n, taken) {
    ok <- rep(TRUE, nrow(cand))
    for (t in seq_len(nrow(taken))) ok <- ok & far(cand, taken[t, ])
    pick <- integer(0)
    for (i in seq_len(n)) {
      j <- which(ok)[which.max(score[ok])]
      if (!length(j)) break
      pick <- c(pick, j)
      ok <- ok & far(cand, cand[j, ])
    }
    pick
  }
  pEI <- greedy(ei, nEI, cand[integer(0), , drop = FALSE])
  pLCB <- greedy(-lcb, nLCB, cand[pEI, , drop = FALSE])
  topi <- ord[seq_len(min(2L * nRep, length(ord)))]
  topi <- topi[order(mult[topi])][seq_len(min(nRep, length(topi)))]
  list(U = rbind(cand[pEI, , drop = FALSE], cand[pLCB, , drop = FALSE], X0[topi, , drop = FALSE]),
       kind = c(rep("ei", length(pEI)), rep("lcb", length(pLCB)), rep("rep", length(topi))),
       bestZ = fstar)
}

#' The objective's arguments, from the fit's `runDEoptim()` arguments
#'
#' What `runDEoptim()` hands `fireSenseUtils::.objfunSpreadFit()`: its arguments that the objective also has, plus
#' `doAssertions` (its `doObjFunAssertions`) and `fitYearSpreadSD` (whether `lower` has `yearSpreadSD`), and the
#' buffer edges it adds to `fireBufferedListDT` before the fit (`fireSenseUtils::addBufferEdge()`). `thresh` is
#' `Inf`: every evaluation is a full one, as in `runDEoptim()`'s re-score, so the GP sees values, not the early
#' stop's fail value.
#'
#' @param deArgs named list, the arguments `fitSpread()` gives `runDEoptim()`.
#' @return a named list of arguments to `.objfunSpreadFit()` (all but `par`).
emulatorObjectiveArgs <- function(deArgs) {
  fn <- fireSenseUtils::.objfunSpreadFit
  args <- deArgs[intersect(names(deArgs), setdiff(names(formals(fn)), "..."))]
  args$fireBufferedListDT <- fireSenseUtils::addBufferEdge(deArgs$fireBufferedListDT, deArgs$landscape)
  args$doAssertions <- deArgs$doObjFunAssertions
  args$fitYearSpreadSD <- "yearSpreadSD" %in% names(deArgs$lower)
  args$thresh <- Inf
  args$verbose <- 0
  args
}

## The objective's data: shipped to the workers once, as runDEoptim() ships them (clusterSetup(objsNeeded)), not
## sent with each evaluation
emulatorDataArgs <- c("landscape", "annualDTx1000", "nonAnnualDTx1000", "fireBufferedListDT", "historicalFires")

#' Run the emulator on the fit's objective, on a cluster built as the fit's is
#'
#' With `cores`, the evaluations run on a cluster from `clusters::clusterSetup()`, asked for `nCoresNeeded` workers
#' on `cores` (booked in clusters' reservations, so other fits see it), with the objective's data shipped once; the
#' cluster is stopped and its booking released when this function returns, as `runDEoptim()`'s is. Without `cores`
#' (`NULL` or `NA`) the evaluations run here, one at a time.
#'
#' @param objArgs `emulatorObjectiveArgs()`.
#' @param lower,upper the bounds, named.
#' @param passingPars the threshold calibration's passing sets, or `NULL`; they centre the box (`emulatorBox()`).
#' @param NP integer; the starting population's size.
#' @param control list; see `emulatorDefaults()`.
#' @param seed integer.
#' @param cores,nCoresNeeded,libPath,logPath,runName as for the fit (`runDEoptim()`).
#' @return the `emulateObjective()` list.
emulatorStartPopulation <- function(objArgs, lower, upper, passingPars, NP, control, seed,
                                    cores = NULL, nCoresNeeded = NP, libPath = .libPaths()[1], logPath = tempdir(),
                                    runName = "") {
  ctl <- emulatorControlOf(control)
  box <- emulatorBox(passingPars, lower, upper, widen = ctl$boxWiden, minWidthFrac = ctl$boxMinWidthFrac)
  if (isTRUE(is.na(cores))) cores <- NULL
  cl <- NULL
  if (!is.null(cores)) {
    withr::local_options(clusters.workerMemoryGB =
                           getOption("clusters.workerMemoryGB", fireSenseUtils::spreadFitWorkerMemoryGB))
    shipped <- list2env(objArgs[emulatorDataArgs])
    built <- clusters::clusterSetup(messagePrefix = paste0("emulator_", runName, "_"), cores = cores,
                                    nCoresNeeded = nCoresNeeded, libPath = libPath[1], logPath = logPath,
                                    objsNeeded = as.list(emulatorDataArgs),
                                    pkgsNeeded = c("data.table", "fireSenseUtils"), envir = shipped,
                                    runName = paste0("emulator_", runName))
    cl <- built$cluster
    if (is.null(ctl$batchSize)) ctl$batchSize <- length(cl)
  }
  message("fireSense_spreadFit: emulator search for the starting population (", length(lower), " parameters, ",
          if (is.null(cl)) "evaluated here" else paste(length(cl), "workers"), ", box from ",
          if (is.null(passingPars) || !NROW(passingPars)) "lower/upper" else paste(NROW(passingPars), "passing sets"), ")")
  emulateObjective(emulatorEvaluator(objArgs, cl), lower, upper, box, NP, ctl, seed)
}

## Evaluates parameter sets (rows) once each with fireSenseUtils::rescorePopulation(), on `cl` (whose workers hold
## the data) or here (the data go with the arguments)
emulatorEvaluator <- function(objArgs, cl = NULL) {
  fnArgs <- if (is.null(cl)) objArgs else objArgs[setdiff(names(objArgs), emulatorDataArgs)]
  function(P, seed) {
    if (!is.null(cl)) cl <- clusters::currentCluster(cl)
    fireSenseUtils::rescorePopulation(P, fn = fireSenseUtils::.objfunSpreadFit, reps = 1L, cl = cl,
                                      seed = seed, fnArgs = fnArgs)$value
  }
}

#' The DEoptim stopping rule of a fit, seeded or not
#'
#' A fit that starts from the emulator's points starts near where an unseeded fit ends, so it stops sooner:
#' `minGenerations` and `window` default to `seededConvergence` (100 each) instead of clusters' defaults. Only those
#' two; `se` and `sdChange` are clusters' as always. Elements of the `convergence` parameter win.
#'
#' @param convergence the `convergence` parameter (a list, or `NULL`).
#' @param seeded logical; the start came from the emulator.
#' @param seededConvergence list; the seeded fit's `window` and `minGenerations`.
#' @return the full rule, `clusters::deoptimConvergence()`.
convergenceOf <- function(convergence, seeded, seededConvergence = list(window = 100L, minGenerations = 100L)) {
  if (seeded) convergence <- c(convergence, seededConvergence[setdiff(names(seededConvergence), names(convergence))])
  clusters::deoptimConvergence(convergence)
}

## The emulator's code, for the Cache key of its result: Cache() digests the arguments, not the functions they reach
emulatorBodies <- function() {
  fns <- c("emulatorDefaults", "emulatorControlOf", "emulatorBox", "emulateObjective", "emulatorScale",
           "fitEmulatorGP", "emulatorConverged", "selectEmulatorBatch", "emulatorEvaluator")
  stats::setNames(lapply(fns, function(f) deparse(get(f, mode = "function"))), fns)
}

#' Whether this fit starts from the emulator
#'
#' `initialPopMethod = "emulator"`, unless the `initialpop` parameter is set (it wins) or there are more parameters
#' than `emulatorMaxParams` (a GP in that many dimensions needs more evaluations than the search saves); each
#' fallback says so and the fit starts as with `"calibration"`.
#'
#' @param Par the module's parameters, `P(sim)`.
#' @param nPar integer; the number of fitted parameters, `length(lower)`.
#' @return logical.
useEmulator <- function(Par, nPar) {
  if (!Par$initialPopMethod %in% c("calibration", "emulator"))
    stop("fireSense_spreadFit: initialPopMethod must be 'calibration' or 'emulator'; got ",
         dQuote(Par$initialPopMethod, FALSE), call. = FALSE)
  if (Par$initialPopMethod != "emulator") return(FALSE)
  if (!is.null(Par$initialpop)) {
    message("fireSense_spreadFit: initialpop is set, so the emulator is not run")
    return(FALSE)
  }
  if (nPar > Par$emulatorMaxParams) {
    message("fireSense_spreadFit: ", nPar, " parameters is more than emulatorMaxParams (", Par$emulatorMaxParams,
            "); starting as with initialPopMethod = 'calibration'")
    return(FALSE)
  }
  TRUE
}
