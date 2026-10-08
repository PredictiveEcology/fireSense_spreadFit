#' Draw parameter sets from the logistic's active range, by rejection
#'
#' Draws sets uniformly between `lower` and `upper` and keeps those `gates` accepts, until `n` are kept
#' or `maxDraws` have been drawn. Uniform draws over the wide bounds mostly saturate the logistic
#' (spreadProb at its asymptote or floor everywhere), which the objective refuses ("Too burny a
#' landscape", "Not spread out enough"); for ELF 14.3 fold 1 only about 7% of them pass. A trial the
#' objective refuses has no usable SNLL, so calibrating on such draws can leave no usable trial and no threshold.
#'
#' The draws are one `runif()` stream from the seed set by the caller, so the result is reproducible; the
#' batch sizes depend only on how many draws were accepted so far, which is deterministic. Drawing
#' stops at the draw that made the `n`th acceptance, so `drawn` counts only the draws made use of.
#'
#' @param n integer; number of parameter sets wanted.
#' @param lower,upper named numeric; the bounds.
#' @param maxDraws integer; cap on the number of draws.
#' @param gates function of a list of parameter sets, returning a logical vector (TRUE = accepted).
#'   Evaluated in forked processes if `cores > 1`.
#' @param cores integer; forks used to evaluate `gates`.
#' @return `list(pars =, drawn =, accepted =)`: the accepted sets (at most `n`), how many draws it took,
#'   and how many were accepted.
drawActivePars <- function(n, lower, upper, maxDraws, gates, cores = 1L) {
  k <- length(lower)
  acc <- list()
  drawn <- 0L
  repeat {
    need <- n - length(acc)
    left <- maxDraws - drawn
    if (need <= 0 || left <= 0) break
    rate <- if (drawn > 0L) max(length(acc) / drawn, 1 / maxDraws) else 1
    size <- as.integer(min(left, max(need, ceiling(1.2 * need / rate))))
    draws <- lapply(seq_len(size), function(i) {
      x <- runif(k, lower, upper)
      names(x) <- names(lower)
      x
    })
    chunks <- split(seq_len(size), rep_len(seq_len(max(1L, min(cores, size))), size))
    ok <- parallel::mclapply(chunks, function(ii) gates(draws[ii]), mc.cores = max(1L, min(cores, length(chunks))))
    if (any(vapply(ok, function(o) !is.logical(o), logical(1))))
      stop("spreadProb gates failed on a draw: ", paste(unique(unlist(lapply(ok, as.character))), collapse = "; "))
    pass <- logical(size)
    for (j in seq_along(chunks)) pass[chunks[[j]]] <- ok[[j]]
    keep <- which(pass)
    if (length(keep) >= need) {
      keep <- keep[seq_len(need)]
      acc <- c(acc, draws[keep])
      drawn <- drawn + keep[need]
      break
    }
    acc <- c(acc, draws[keep])
    drawn <- drawn + size
  }
  list(pars = acc, drawn = drawn, accepted = length(acc))
}
