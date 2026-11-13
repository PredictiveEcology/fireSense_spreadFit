#' The spread objective's gates, as a function of parameter sets
#'
#' Builds the screen `feasibleInitialPopulation()` applies to its draws:
#' `fireSenseUtils::spreadProbGates()` on the covariates of the first block of fire years, which is
#' the test the objective itself applies before simulating any fire. A parameter set that fails it is
#' one the objective refuses (its fail value, 1e6).
#'
#' @param covs `mod$covsX1000`, or the subset of its years being fitted.
#' @param formulaToFit,covMinMax,covCentre,climateRef,mutuallyExclusive,maxFireSpread,escapeSizeHa,landscape,link,powerLawGamma,powerLawShape,lowerAsymptote,richardsNu
#'   as given to `fireSenseUtils::runDEoptim()` for the same fit.
#' @return a function of a list of parameter sets (named numeric vectors, as `drawActivePars()` draws
#'   them) returning a logical vector, TRUE for a set that passes; with `pass = FALSE` the
#'   `spreadProbGates()` data.frame (`pass`, `burny`, `notSpread`, `medianOK`) instead.
spreadGates <- function(covs, formulaToFit, covMinMax, covCentre, climateRef, mutuallyExclusive,
                        maxFireSpread, escapeSizeHa, landscape, link, powerLawGamma, powerLawShape,
                        lowerAsymptote = fireSenseUtils::spreadProbFloor, richardsNu = fireSenseUtils::spreadRichardsNu) {
  force(covs); force(formulaToFit); force(covMinMax); force(covCentre); force(climateRef)
  force(mutuallyExclusive); force(maxFireSpread); force(escapeSizeHa); force(landscape); force(link)
  force(powerLawGamma); force(powerLawShape); force(lowerAsymptote); force(richardsNu)
  function(parList, pass = TRUE) {
    g <- fireSenseUtils::spreadProbGates(
      parList, annualDTx1000 = covs$annualDTx1000, nonAnnualDTx1000 = covs$nonAnnualDTx1000,
      historicalFires = covs$historicalFires, formulaToFit = formulaToFit,
      covMinMax = covMinMax, covCentre = covCentre, climateRef = climateRef,
      mutuallyExclusive = mutuallyExclusive, maxFireSpread = maxFireSpread, escapeSizeHa = escapeSizeHa,
      landscape = landscape, link = link, powerLawGamma = powerLawGamma, powerLawShape = powerLawShape,
      lowerAsymptote = lowerAsymptote, richardsNu = richardsNu, doAssertions = FALSE)
    if (pass) g$pass else g
  }
}

#' A DEoptim starting population the objective can fit from
#'
#' DEoptim's own start is uniform between the bounds, which are deliberately wide (sign only). For an
#' ELF whose feasible region is a sliver of that box, the start has no member the objective accepts, and
#' DE never finds one: a trial is kept when it is no worse than its parent, and every value is the
#' fail value (ELF 6.6.2, 2026-10-20: 0 of 2000 uniform draws passed the gates, with two climate terms
#' each drawn over 0-50). This draws by rejection (`drawActivePars()`) and fills the rest of the population with uniform draws: DE's p-best strategy pulls the
#' population towards its accepted members, so a few are enough to start from.
#'
#' @param n integer; rows wanted. At least the population size the fit may get (`clusters::clusterSetup()`
#'   uses the first `NP` rows, `NP` being the number of workers).
#' @param lower,upper the bounds, named.
#' @param gates a function of a list of parameter sets, see `spreadGates()`.
#' @param seed integer; `withr::with_seed()`, so the population is the same on every call for an ELF.
#' @param minFeasible integer; stop with the gate counts when fewer sets pass after `maxDrawsPerMember * n`
#'   draws. 4 is DEoptim's smallest population.
#' @param maxDrawsPerMember integer; cap on the draws per member wanted.
#' @param cores integer; forks for `gates`.
#' @return an `n` x `length(lower)` matrix, columns named as `lower`, the accepted sets first; attribute
#'   `feasible` is how many rows were accepted.
feasibleInitialPopulation <- function(n, lower, upper, gates, seed, minFeasible = 4L,
                                      maxDrawsPerMember = 100L, cores = 1L) {
  n <- as.integer(n)
  withr::with_seed(seed, {
    active <- drawActivePars(n, lower, upper, maxDraws = as.integer(maxDrawsPerMember * n), gates, cores = cores)
    if (active$accepted < minFeasible) {
      counts <- gates(lapply(seq_len(200L), function(i) stats::setNames(runif(length(lower), lower, upper), names(lower))),
                      pass = FALSE)
      stop("fireSense_spreadFit: only ", active$accepted, " of ", active$drawn, " parameter sets drawn between the bounds ",
           "pass the objective's spreadProb gates (fewer than ", minFeasible, "); there is no population to fit from. ",
           "Of 200 further draws: burny ", sum(counts$burny), ", not spread out enough ", sum(counts$notSpread),
           ", median in range ", sum(counts$medianOK), ". Look at the covariate scales and bounds of this ELF.",
           call. = FALSE)
    }
    pop <- do.call(rbind, c(active$pars, uniformSets(n - active$accepted, lower, upper)))
  })
  colnames(pop) <- names(lower)
  message("Starting population: ", active$accepted, " of ", n, " members pass the spreadProb gates (",
          active$drawn, " draws); the rest are uniform between the bounds")
  attr(pop, "feasible") <- active$accepted
  pop
}

## `m` parameter sets drawn uniformly between the bounds, as a list of named vectors
uniformSets <- function(m, lower, upper)
  lapply(seq_len(m), function(i) stats::setNames(runif(length(lower), lower, upper), names(lower)))

## `DEoptimControl` with a starting population added; unchanged (same cache key) without one
deoptimControlStart <- function(control, initialpop = NULL) {
  if (is.null(initialpop)) control else c(control, list(initialpop = initialpop))
}

#' The population DEoptim starts from
#'
#' The `initialpop` parameter's rows first, then uniform draws up to `n` rows, when it is given (it used to be dropped: `runDEoptim()` ignores its own `initialpop`
#' argument, so a population seeded from elsewhere never reached DEoptim, 2026-10-26); otherwise `NULL`,
#' DEoptim's own draw. The rows were drawn between `lower` and `upper`, so one outside them is a changed bound
#' or a bug and stops, rather than being clamped.
#'
#' @param initialpop the `initialpop` parameter: a matrix with the columns of `lower`, or `NULL`.
#' @param n rows of the population (at least the rows of `initialpop`).
#' @param lower,upper the bounds, named.
#' @param seed seed for the uniform draws.
#' @return A matrix or `NULL`.
startPopulation <- function(initialpop, n, lower, upper, seed) {
  if (is.null(initialpop)) return(NULL)
  pop <- as.matrix(initialpop)
  if (!is.null(colnames(pop)) && !identical(colnames(pop), names(lower)))
    stop("initialpop's columns (", paste(colnames(pop), collapse = ", "), ") are not names(lower) in that order",
         call. = FALSE)
  if (ncol(pop) != length(lower)) stop("initialpop has ", ncol(pop), " columns; lower has ", length(lower), call. = FALSE)
  colnames(pop) <- names(lower)
  if (any(sweep(pop, 2L, lower) < 0) || any(sweep(pop, 2L, upper) > 0))
    stop("fireSense_spreadFit: a row of initialpop is outside lower/upper", call. = FALSE)
  uniform <- do.call(rbind, withr::with_seed(seed, uniformSets(max(n, nrow(pop)), lower, upper)))
  out <- rbind(pop, uniform)[seq_len(nrow(uniform)), , drop = FALSE]
  colnames(out) <- names(lower)
  out
}
