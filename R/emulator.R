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
#'   - `minEvals` (600): the global phase does not stop for lack of improvement before this many evaluations;
#'   - `globalWindow` (6), `relTol` (0.002): the global phase stops when the GP's best mean (at an evaluated point)
#'     improved by less than `relTol` times itself over its last `globalWindow` batches (the prototype: 2 objective
#'     units at ~938);
#'   - `batchSize` (`NULL`): evaluations per batch; `NULL` is the cluster's number of workers, or 48 without one;
#'   - `floorFrac` (0.8), `globalCapRatio` (10), `localCapRatio` (1.18): the GP models
#'     `log(min(y, cap) - floorFrac * min(y))`, with `cap` `globalCapRatio` times the best value seen in the global
#'     phase and `localCapRatio` times the GP's best mean in the local phase, so failures (the objective's fail
#'     value, and runaway sets scored in the millions) sit on a plateau instead of dominating the fit;
#'   - `nCandidates` (4000), `nPerturb` (3000), `perturbSD` (c(0.005, 0.015, 0.04, 0.1)), `nTop` (20): each batch
#'     is chosen from `nCandidates` Latin-hypercube points (in the trust region in the local phase) and `nPerturb`
#'     Gaussian perturbations, at each `perturbSD` (unit-cube scale), of the `nTop` best points by GP mean;
#'   - `exclusion` (0.25): batch members are at least this far apart in lengthscale units;
#'   - `lcbSD` (2): the exploration picks minimise the GP mean minus this many standard deviations;
#'   - `logSpanMin` (20): the search space (`emulatorSpace()`) puts a parameter on a signed log1p scale when one of
#'     its bounds is 0 and the other is at least this far from 0; the others are linear;
#'   - the local phase, trust-region rounds of `trustRoundBatches` (2) batches each: the region is the box around
#'     the best `localTopFrac` (0.1) of all distinct evaluated points by GP mean, widened by `localWiden` (0.5) of
#'     its width each side, times the round's scale (`emulatorTrustBox()`). After each round (`emulatorTrustStep()`),
#'     an improvement of the GP's best mean by more than `trustTolSD` (1) times the objective's noise SD
#'     (`emulatorNoiseSD()`), and by more than `relTol` times the best mean (the floor while there are few or no
#'     replicates near it), is a success, else
#'     a failure; `trustSuccessTol` (1) successes in a row multiply the scale by `trustShrink` (0.5),
#'     a failure multiplies it by `trustExpand` (1.25, at most 1), and `trustFailureTol` (2) failures in a row stop
#'     the search. It also stops when the region's widest side is narrower than `trustMinWidth` (0.01 of the unit
#'     cube; not relative to the GP's lengthscales, which reach their upper bound when the local cap flattens the
#'     response), or at `maxEvals`;
#'   - `popPoolFrac` (0.2), `popMinDist` (1): the starting population is picked from the best `popPoolFrac` of the
#'     distinct evaluated points (by their mean value; at least `NP` of them), best first, skipping a point closer
#'     than `popMinDist` GP lengthscales to one already picked, the distance halved until `NP` are picked
#'     (`emulatorSpreadPopulation()`), so separate good regions are all in it;
#'   - `thetaBounds` (c(0.01, 5)), `gpMaxit` (200), `gpThreads` (6): the GP's lengthscale bounds (unit cube),
#'     optimiser iterations (`hetGP::mleHetGP()`) and BLAS threads (the rest of the search runs on one).
emulatorDefaults <- function() {
  list(designPerParam = 20L, designRepPoints = 4L, designReps = 6L,
       globalMaxEvals = 1000L, maxEvals = 1600L, minEvals = 600L,
       globalWindow = 6L, relTol = 0.002, batchSize = NULL,
       floorFrac = 0.8, globalCapRatio = 10, localCapRatio = 1.18,
       nCandidates = 4000L, nPerturb = 3000L, perturbSD = c(0.005, 0.015, 0.04, 0.1), nTop = 20L,
       exclusion = 0.25, lcbSD = 2, logSpanMin = 20,
       trustRoundBatches = 2L, localTopFrac = 0.1, localWiden = 0.5, trustTolSD = 1,
       trustShrink = 0.5, trustExpand = 1.25, trustSuccessTol = 1L, trustFailureTol = 2L, trustMinWidth = 0.01,
       popPoolFrac = 0.2, popMinDist = 1,
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

#' The space the emulator searches: the full bounds, as a unit cube
#'
#' Every parameter spans `lower` to `upper` (the search is not restricted). A parameter whose interval has 0 as one
#' bound and whose other bound is at least `logSpanMin` from 0 (the fuel coefficients, `0`-`100`, and `youngAge`,
#' `-50`-`0`; also `CMD`, `0`-`50`) is searched on a signed log1p scale, `sign(x) * log1p(|x|)`, so the values
#' near its 0 bound, a sign constraint where fitted values sit (1-2 against a bound of 100), get as much of the
#' cube as the far end. The others (bounds off 0, or a narrow range such as `yearSpreadSD`) are linear. The
#' covariates are scaled, so a coefficient's units are comparable across parameters and `logSpanMin` is in them.
#'
#' @param lower,upper the bounds, named.
#' @param logSpanMin numeric; see above.
#' @return a list: `logPars` (names), `toUnit(P)` and `fromUnit(U)` (matrices, columns `names(lower)`;
#'   `fromUnit()` stays within `lower`/`upper`).
emulatorSpace <- function(lower, upper, logSpanMin = 20) {
  parNames <- names(lower)
  isLog <- (lower == 0 | upper == 0) & pmax(abs(lower), abs(upper)) >= logSpanMin
  fwd <- function(P) { P[, isLog] <- sign(P[, isLog]) * log1p(abs(P[, isLog])); P }
  inv <- function(Z) { Z[, isLog] <- sign(Z[, isLog]) * expm1(abs(Z[, isLog])); Z }
  tlo <- fwd(rbind(lower))[1L, ]; twidth <- fwd(rbind(upper))[1L, ] - tlo
  list(logPars = parNames[isLog],
       toUnit = function(P) sweep(sweep(fwd(as.matrix(P)[, parNames, drop = FALSE]), 2L, tlo), 2L, twidth, "/"),
       fromUnit = function(U) {
         P <- inv(sweep(sweep(as.matrix(U), 2L, twidth, "*"), 2L, tlo, "+"))
         P <- sweep(sweep(P, 2L, lower, pmax), 2L, upper, pmin) # expm1(log1p(x)) can be off by one ulp
         colnames(P) <- parNames
         P
       })
}

#' The local phase's trust region
#'
#' The box around the best `topFrac` of the points `X` by `score` (at least 2), widened by `widen` times its width
#' each side, then scaled by `scale` about its centre, and clipped to the unit cube.
#'
#' @param X matrix; distinct points, unit cube.
#' @param score numeric; one per row of `X`, lower is better (the GP mean).
#' @param topFrac,widen,scale numeric.
#' @return a 2-row matrix (`lower`, `upper`).
emulatorTrustBox <- function(X, score, topFrac, widen, scale = 1) {
  top <- X[order(score)[seq_len(min(nrow(X), max(2L, ceiling(topFrac * nrow(X)))))], , drop = FALSE]
  lo <- apply(top, 2L, min); hi <- apply(top, 2L, max)
  mid <- (lo + hi) / 2; half <- (hi - lo) * (0.5 + widen) * scale
  rbind(lower = pmax(mid - half, 0), upper = pmin(mid + half, 1))
}

#' One trust-region round's outcome
#'
#' @param state list: `scale`, `successes`, `failures`.
#' @param improvement numeric; the GP's best mean before the round minus after it.
#' @param tol numeric; an improvement above it is a success.
#' @param ctl `emulatorControlOf()`.
#' @return `state` updated, with `stop` (logical) and `reason`.
emulatorTrustStep <- function(state, improvement, tol, ctl) {
  state$stop <- FALSE; state$reason <- NA_character_
  if (isTRUE(improvement > tol)) {
    state$failures <- 0L
    state$successes <- state$successes + 1L
    if (state$successes >= ctl$trustSuccessTol) {
      state$scale <- state$scale * ctl$trustShrink
      state$successes <- 0L
    }
  } else {
    state$successes <- 0L
    state$failures <- state$failures + 1L
    if (state$failures >= ctl$trustFailureTol) {
      state$stop <- TRUE; state$reason <- "improvement below noise"
    } else state$scale <- min(1, state$scale * ctl$trustExpand)
  }
  state
}

#' The objective's noise near the best points
#'
#' The pooled within-point SD of the replicated evaluations among the best `topFrac` of the distinct points
#' (by their mean value), or `NA` without replicates there.
#'
#' @param U matrix; every evaluated point (unit cube), one row per evaluation.
#' @param y numeric; the evaluations.
#' @param topFrac numeric.
#' @return numeric.
emulatorNoiseSD <- function(U, y, topFrac) {
  ok <- !is.na(y)
  keys <- emulatorPointKeys(U)[ok]; y <- y[ok]
  m <- tapply(y, keys, mean); v <- tapply(y, keys, stats::var); n <- tapply(y, keys, length)
  top <- order(m)[seq_len(min(length(m), max(2L, ceiling(topFrac * length(m)))))]
  top <- top[n[top] > 1L]
  if (!length(top)) return(NA_real_)
  sqrt(sum((n[top] - 1) * v[top]) / sum(n[top] - 1))
}

## One string per row of `U`, the same for the replicates of a point
emulatorPointKeys <- function(U) apply(U, 1L, function(r) paste(sprintf("%a", r), collapse = ","))

#' Search a noisy objective with a Gaussian-process emulator and return its best evaluated points
#'
#' On the full bounds, in `emulatorSpace()`: a space-filling design (`lhs::maximinLHS()`, plus replicates), then
#' batches chosen on a GP (`hetGP::mleHetGP()`, Matern 5/2) fitted to every evaluation so far, on a capped log
#' scale (see `emulatorDefaults()`): expected improvement over the GP's best mean (greedy, with an exclusion
#' distance), a lower confidence bound for exploration, and replicates of the best points. A global phase, until
#' the GP's best mean stops improving, then trust-region rounds (`emulatorTrustBox()`, `emulatorTrustStep()`),
#' each re-centred on the best of all the points so far, shrinking while they improve by more than the
#' objective's noise, until they stop improving or the budget is spent.
#'
#' @param evaluate a function of a matrix of parameter sets (one per row, columns `names(lower)`) and an integer
#'   `seed`, returning one objective value per row (`emulatorEvaluator()`).
#' @param lower,upper the bounds, named.
#' @param NP integer; how many points to return.
#' @param control list; see `emulatorDefaults()`.
#' @param seed integer; the design, the candidates and the evaluations' seeds follow from it.
#' @return a list: `pop`, `NP` distinct evaluated points, good and spread out (`emulatorSpreadPopulation()`; fewer
#'   if fewer were evaluated), columns `names(lower)`; `objFunVal`, the means of their evaluations; `evals`, every
#'   evaluation (`batch`, `kind`, `value` and the parameters); `status`, one row per batch (with the trust `round`, its `scale` and its
#'   box's widest side, `boxWidth`); `stopReason`.
emulateObjective <- function(evaluate, lower, upper, NP, control = list(), seed = 1L) {
  ctl <- emulatorControlOf(control)
  d <- length(lower)
  space <- emulatorSpace(lower, upper, ctl$logSpanMin)
  B <- if (is.null(ctl$batchSize)) 48L else as.integer(ctl$batchSize)
  localSingleThreadBLAS() # the evaluations, and the GP's predictions; its fits take `gpThreads`
  withr::with_seed(seed, {
    U <- lhs::maximinLHS(as.integer(ctl$designPerParam * d), d)
    repRows <- unique(round(seq(1, nrow(U), length.out = ctl$designRepPoints + 1L)))[seq_len(ctl$designRepPoints)]
    U <- rbind(U, U[rep(repRows, each = ctl$designReps), , drop = FALSE])
    kind <- c(rep("design", nrow(U) - length(repRows) * ctl$designReps), rep("rep", length(repRows) * ctl$designReps))
    y <- evaluate(space$fromUnit(U), seed = seed)
    batch <- rep(0L, nrow(U))
    status <- NULL
    k <- 0L
    theta <- 1 # the latest GP's lengthscales
    ## the GP on every evaluation so far; `local`: capped near the best (see emulatorDefaults())
    fitGP <- function(local) {
      best <- if (is.null(status)) min(y, na.rm = TRUE) else min(status$bestMean)
      cap <- if (local) ctl$localCapRatio * best else ctl$globalCapRatio * min(y, na.rm = TRUE)
      sc <- emulatorScale(y, cap, ctl$floorFrac)
      list(m = fitEmulatorGP(U, sc$z, ctl), sc = sc)
    }
    ## evaluate a batch and record it; trustScale NULL is the global phase
    evaluateBatch <- function(sel, sc, trustScale = NULL, round = NA_integer_) {
      k <<- k + 1L
      local <- !is.null(trustScale)
      theta <<- sel$theta
      new <- evaluate(space$fromUnit(sel$U), seed = seed + 1e6 * k)
      U <<- rbind(U, sel$U); y <<- c(y, new); kind <<- c(kind, sel$kind); batch <<- c(batch, rep(k, nrow(sel$U)))
      status <<- rbind(status, data.table(batch = k, phase = if (local) "local" else "global", round = round,
                                          scale = if (local) trustScale else NA_real_,
                                          boxWidth = if (local) max(sel$box[2L, ] - sel$box[1L, ]) else NA_real_,
                                          nEval = nrow(U), bestMean = sc$untr(sel$bestZ),
                                          batchBest = min(new, na.rm = TRUE)))
      message("emulator batch ", k, if (local) paste0(" (trust round ", round, ", scale ", signif(trustScale, 3), ")")
              else " (global)", ": ", nrow(U), " evaluations; GP best mean ",
              signif(utils::tail(status$bestMean, 1L), 6), "; best single value ", signif(min(y, na.rm = TRUE), 6))
    }
    while (nrow(U) < min(ctl$globalMaxEvals, ctl$maxEvals) && !emulatorConverged(status, nrow(U), ctl)) {
      g <- fitGP(local = FALSE)
      evaluateBatch(selectEmulatorBatch(g$m, U, B, NULL, ctl), g$sc)
    }
    ## trust-region rounds. A round is judged by the GP fitted after its evaluations: its best mean against the
    ## one at the round's start (both local-phase GPs, on the same capped scale)
    trust <- list(scale = 1, successes = 0L, failures = 0L, stop = FALSE, reason = "maxEvals")
    round <- 1L; inRound <- 0L; before <- NULL
    while (nrow(U) < ctl$maxEvals) {
      g <- fitGP(local = TRUE)
      sel <- selectEmulatorBatch(g$m, U, B, trust$scale, ctl)
      now <- g$sc$untr(sel$bestZ)
      if (is.null(before)) before <- now
      if (inRound >= ctl$trustRoundBatches) {
        tol <- max(ctl$trustTolSD * emulatorNoiseSD(U, y, ctl$localTopFrac), ctl$relTol * now, na.rm = TRUE)
        trust <- emulatorTrustStep(trust, before - now, tol, ctl)
        message("emulator trust round ", round, ": the GP's best mean improved by ", signif(before - now, 3),
                " (tolerance ", signif(tol, 3), ")", if (!trust$stop) paste0("; next scale ", signif(trust$scale, 3)))
        if (trust$stop) break
        round <- round + 1L; inRound <- 0L; before <- now
        sel <- selectEmulatorBatch(g$m, U, B, trust$scale, ctl)
        if (max(sel$box[2L, ] - sel$box[1L, ]) < ctl$trustMinWidth) {
          trust$reason <- "trust region narrower than trustMinWidth"
          break
        }
      }
      evaluateBatch(sel, g$sc, trust$scale, round)
      inRound <- inRound + 1L
    }
  })
  ok <- which(!is.na(y))
  keys <- emulatorPointKeys(U)[ok]
  rows <- ok[!duplicated(keys)]                                           # each distinct point once
  means <- as.numeric(tapply(y[ok], factor(keys, levels = unique(keys)), mean))
  pick <- emulatorSpreadPopulation(U[rows, , drop = FALSE], means, theta, NP, ctl$popPoolFrac, ctl$popMinDist)
  list(pop = space$fromUnit(U[rows[pick], , drop = FALSE]), objFunVal = means[pick],
       evals = data.table(batch = batch, kind = kind, value = y, space$fromUnit(U)), status = status,
       stopReason = trust$reason)
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

## The GP fitted to every evaluation, on the unit cube (hetGP pools replicates)
fitEmulatorGP <- function(U, z, ctl) {
  localSingleThreadBLAS(threads = ctl$gpThreads)
  d <- ncol(U)
  hetGP::mleHetGP(X = U, Z = z, lower = rep(ctl$thetaBounds[1L], d), upper = rep(ctl$thetaBounds[2L], d),
                  covtype = "Matern5_2", maxit = ctl$gpMaxit, settings = list(trace = -1))
}

## Whether the global phase has stopped improving: at least `minEvals` evaluations, more than `globalWindow`
## batches, and the GP's best mean improved by less than `relTol` times itself over the last `globalWindow` of them
emulatorConverged <- function(status, nEval, ctl) {
  if (is.null(status) || nEval < ctl$minEvals) return(FALSE)
  b <- status$bestMean
  n <- length(b)
  n > ctl$globalWindow && b[n - ctl$globalWindow] - b[n] < ctl$relTol * b[n]
}

#' One batch of points to evaluate, chosen on the GP
#'
#' Expected improvement over the GP's best mean at an evaluated point (greedy, each pick at least
#' `exclusion` lengthscales from the others), then a lower confidence bound (exploration), then replicates of the
#' best evaluated points (fewest replicates first). In the local phase the Latin-hypercube candidates are drawn in
#' the trust region (`emulatorTrustBox()`) and more replicates are taken.
#'
#' @param m the GP (`fitEmulatorGP()`).
#' @param U every evaluated point so far, unit cube.
#' @param B integer; batch size.
#' @param trustScale numeric; the trust region's scale in the local phase, `NULL` in the global phase.
#' @param ctl `emulatorControlOf()`.
#' @return list: `U` (the batch, unit cube), `kind` ("ei", "lcb", "rep"), `bestZ` (the GP's best mean at an
#'   evaluated point, response scale), `box` (the trust region, or `NULL`), `theta` (the GP's lengthscales).
selectEmulatorBatch <- function(m, U, B, trustScale, ctl) {
  local <- !is.null(trustScale)
  nRep <- round(B / if (local) 4 else 6); nLCB <- round(B / 8); nEI <- B - nRep - nLCB
  d <- ncol(U)
  keys <- emulatorPointKeys(U)
  X0 <- U[!duplicated(keys), , drop = FALSE]
  mult <- tabulate(match(keys, keys[!duplicated(keys)]), nrow(X0))
  pm <- predict(m, X0)$mean
  ord <- order(pm)
  top <- X0[ord[seq_len(min(ctl$nTop, nrow(X0)))], , drop = FALSE]
  cand <- lhs::randomLHS(ctl$nCandidates, d)
  box <- NULL
  if (local) {
    box <- emulatorTrustBox(X0, pm, ctl$localTopFrac, ctl$localWiden, trustScale)
    cand <- sweep(sweep(cand, 2L, box[2L, ] - box[1L, ], "*"), 2L, box[1L, ], "+")
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
       bestZ = fstar, box = box, theta = theta)
}

#' The starting population: good evaluated points, spread out
#'
#' From the best `poolFrac` of the distinct points by `value` (at least `NP`), best first, each point that is at
#' least `minDist` from every point already picked, the distance in GP lengthscales (`theta`) of the unit cube.
#' While fewer than `NP` are picked, the distance is halved and the pool is gone through again for more, keeping
#' the points already picked (so each good region gets a point before any gets a second; at 0 every pool point
#' qualifies).
#'
#' @param X matrix; distinct points, unit cube.
#' @param value numeric; one per row of `X` (the mean of its evaluations), lower is better.
#' @param theta numeric; the GP's lengthscales, one per column of `X` (or one for all).
#' @param NP integer.
#' @param poolFrac,minDist numeric.
#' @return integer; rows of `X`, in the order picked, `min(NP, nrow(X))` of them.
emulatorSpreadPopulation <- function(X, value, theta, NP, poolFrac, minDist) {
  theta <- rep_len(theta, ncol(X))
  pool <- order(value)[seq_len(min(length(value), max(NP, ceiling(poolFrac * length(value)))))]
  n <- min(NP, length(pool))
  pick <- pool[1L]
  dist <- minDist
  while (length(pick) < n) {
    for (i in setdiff(pool, pick)) {
      if (length(pick) >= n) break
      if (all(sqrt(colSums(((t(X[pick, , drop = FALSE]) - X[i, ]) / theta)^2)) >= dist)) pick <- c(pick, i)
    }
    dist <- if (dist > minDist / 1024) dist / 2 else 0
  }
  pick
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
#' @param NP integer; the starting population's size.
#' @param control list; see `emulatorDefaults()`.
#' @param seed integer.
#' @param cores,nCoresNeeded,libPath,logPath,runName as for the fit (`runDEoptim()`).
#' @return the `emulateObjective()` list.
emulatorStartPopulation <- function(objArgs, lower, upper, NP, control, seed,
                                    cores = NULL, nCoresNeeded = NP, libPath = .libPaths()[1], logPath = tempdir(),
                                    runName = "") {
  ctl <- emulatorControlOf(control)
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
          if (is.null(cl)) "evaluated here" else paste(length(cl), "workers"), ", between lower and upper; on a log1p ",
          "scale: ", paste(emulatorSpace(lower, upper, ctl$logSpanMin)$logPars, collapse = ", "), ")")
  emulateObjective(emulatorEvaluator(objArgs, cl), lower, upper, NP, ctl, seed)
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
#' `minGenerations` and `window` default to `seededConvergence` (200 each) instead of clusters' defaults. Only those
#' two; `se` and `sdChange` are clusters' as always. Elements of the `convergence` parameter win.
#'
#' @param convergence the `convergence` parameter (a list, or `NULL`).
#' @param seeded logical; the start came from the emulator.
#' @param seededConvergence list; the seeded fit's `window` and `minGenerations`.
#' @return the full rule, `clusters::deoptimConvergence()`.
convergenceOf <- function(convergence, seeded, seededConvergence = list(window = 200L, minGenerations = 200L)) {
  if (seeded) convergence <- c(convergence, seededConvergence[setdiff(names(seededConvergence), names(convergence))])
  clusters::deoptimConvergence(convergence)
}

## The emulator's code, for the Cache key of its result: Cache() digests the arguments, not the functions they reach
emulatorBodies <- function() {
  fns <- c("emulatorDefaults", "emulatorControlOf", "emulatorSpace", "emulatorTrustBox", "emulatorTrustStep",
           "emulatorNoiseSD", "emulatorPointKeys", "emulateObjective", "emulatorScale", "fitEmulatorGP",
           "emulatorConverged", "selectEmulatorBatch", "emulatorSpreadPopulation", "emulatorEvaluator")
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
