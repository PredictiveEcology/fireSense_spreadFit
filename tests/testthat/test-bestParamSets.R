## The ledger's parameter sets for an ELF.
##
## FireSense, 2026-09-17: the ledger's "5 best" were five copies of one member. The run event took the 5
## GENERATIONS with the lowest best value, and those generations all hold the same frozen best member --
## itself usually a lucky draw (ranked 1st to 6th of 60 by replicated mean in eight fits). The ledger also
## never stored covMinMax_spread, which prediction needs to rescale covariates.

## a DEoptim result: one element per generation; the last holds the final population and its values
fakeDE <- function(pop, popval, rescore = NULL) {
  gen <- list(member = list(pop = pop, popval = popval, bestvalit = min(popval)))
  DE <- rep(list(gen), 3)
  if (!is.null(rescore)) attr(DE, "finalRescore") <- rescore
  DE
}
pop <- cbind(c(1, 2, 3, 4, 5, 6, 7), 0)
parNames <- c("p1", "p2")

test_that("with a re-score, the n distinct members with the lowest replicated means are recorded", {
  skip_if_not("bestByReplicatedMean" %in% getNamespaceExports("fireSenseUtils"))
  ## DEoptim holds member 1 as best (a lucky 0.1); replicated means rank members 3, 2, 4, 5, 6 first
  rescore <- data.table::data.table(member = rep(1:7, each = 2), rep = rep(1:2, 7),
                                    value = rep(c(9, 2, 1, 3, 4, 5, 8), each = 2))
  out <- bestParamSets(fakeDE(pop, c(0.1, 5, 6, 7, 8, 9, 10), rescore), parNames)
  expect_equal(out$params$p1, c(3, 2, 4, 5, 6))
  expect_equal(out$objFunVal, c(1, 2, 3, 4, 5))
  expect_identical(names(out$params), parNames)
})

test_that("without a re-score, the final population's own values are used, one row per distinct member", {
  skip_if_not("bestByReplicatedMean" %in% getNamespaceExports("fireSenseUtils"))
  dupPop <- rbind(pop, pop[1, ])             # member 1 twice
  out <- bestParamSets(fakeDE(dupPop, c(1, 2, 3, 4, 5, 6, 7, 1)), parNames)
  expect_identical(nrow(unique(out$params)), 5L)
  expect_equal(out$params$p1, c(1, 2, 3, 4, 5))
})

test_that("the run event records bestParamSets() and covMinMax_spread in the ledger row", {
  mainFiles <- list.files(moduleRoot, pattern = "\\.R$", full.names = TRUE)
  mainFile <- mainFiles[vapply(mainFiles, function(f)
    any(grepl("defineModule\\(", readLines(f, warn = FALSE))), logical(1))]
  src <- paste(readLines(mainFile, warn = FALSE), collapse = "\n")
  expect_match(src, "bestParamSets(DE,", fixed = TRUE)
  rowSrc <- paste(readLines(file.path(moduleRoot, "R", "spreadFitLedgerRow.R"), warn = FALSE), collapse = "\n")
  expect_match(rowSrc, "I(list(sim$covMinMax_spread))", fixed = TRUE)   # the builder `run` and the folds share
  expect_false(grepl("head(sim$DE, 5)", src, fixed = TRUE))
})
