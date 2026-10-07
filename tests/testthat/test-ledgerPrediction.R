## fireSense_spreadPredict reads a ledger row's parameters by splitting `sa$params[[1]][ind, ]` into
## `covPars` (names matching the formula's covariates) and `logisticPars` (everything else), then
## calls `fireSenseUtils::logisticAll(logisticPars, mat, covPars, lowerSpreadProb)`
## (fireSense_spreadPredict.R). That split is by NAME, so an OLD ledger row (fitted before hillSlope1
## was fixed, carrying its own fitted `hillSlope1`) and a NEW one (via `addFixedParsToLedger()`,
## `hillSlope1 = 1`, `inflectionPoint1 = 1`) both predict correctly, each with its own value -- reproduced here without a
## full SpreadPredict simList.

predictRow <- function(par) {
  mat <- matrix(c(0.1, 0.4, 0.7, 1), ncol = 1, dimnames = list(NULL, "cov"))
  covPars <- par[intersect(names(par), "cov")]
  logisticPars <- par[setdiff(names(par), names(covPars))]
  fireSenseUtils::logisticAll(logisticPars, mat, covPars, lowerSpreadProb = 0.13)
}

test_that("an old ledger row keeps predicting with its own fitted hillSlope1", {
  oldRow <- c(maxAsymptote = 0.27, hillSlope1 = 0.5, inflectionPoint1 = 4, cov = 2)
  atFixed <- c(maxAsymptote = 0.27, hillSlope1 = 1, inflectionPoint1 = 4, cov = 2)
  expect_false(isTRUE(all.equal(predictRow(oldRow), predictRow(atFixed))))
  ## matches calling the link directly with hillSlope1 = 0.5
  expect_identical(predictRow(oldRow),
                   fireSenseUtils::logistic3p(matrix(c(0.1, 0.4, 0.7, 1), ncol = 1) %*% 2,
                                              c(0.27, 0.5, 4), par1 = 0.13))
})

test_that("a new ledger row (addFixedParsToLedger()) predicts with hillSlope1 = 1 and inflectionPoint1 = 1", {
  paramsBest <- data.table::data.table(maxAsymptote = 0.27, cov = 2)   # as bestParamSets() returns it
  ## as fireSense_spreadPredict.R reads a ledger row: `sa$params[[1]][ind,] |> as.vector() |> unlist()`
  newRow <- addFixedParsToLedger(paramsBest)[1, ] |> as.vector() |> unlist()
  expect_identical(names(newRow), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "cov"))
  expect_equal(unname(newRow[c("hillSlope1", "inflectionPoint1")]), c(1, 1))
  expect_identical(predictRow(newRow),
                   fireSenseUtils::logistic3p(matrix(c(0.1, 0.4, 0.7, 1), ncol = 1) %*% 2,
                                              c(0.27, 1, 1), par1 = 0.13))
})

test_that("the ledger row from a short par has both fixed values in position, also with upperTail1", {
  short <- data.table::data.table(maxAsymptote = c(0.27, 0.26), upperTail1 = c(-0.3, -0.2), cov = c(2, 3))
  row <- addFixedParsToLedger(short)
  expect_identical(names(row), c("maxAsymptote", "hillSlope1", "inflectionPoint1", "upperTail1", "cov"))
  expect_identical(row$hillSlope1, c(1, 1))
  expect_identical(row$inflectionPoint1, c(1, 1))
  expect_identical(row$cov, c(2, 3))
  expect_identical(names(short), c("maxAsymptote", "upperTail1", "cov"))   # the input is not modified
})
