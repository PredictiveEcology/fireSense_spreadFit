## `skipUnusedFitPrep`: a job that only predicts from the ledger's fit does not schedule `spreadFitPrepare`
## (it makes covariates only the fit, `crossValidate` and the held-out fold read). Toy inputs; `init` only.

ledgerObjs <- function() list(studyAreaWithSpreadParams = toyLedger(c("9.8", "9.9")))
skipInit <- function(params = list(), objects = ledgerObjs())
  queued(runEvents(toySim(c(list(skipUnusedFitPrep = TRUE), params), objects), "init"))

test_that("by default, a predict-only job still prepares the covariates (unchanged)", {
  expect_identical(queued(runEvents(toySim(list(doFit = FALSE), ledgerObjs()), "init")), "spreadFitPrepare")
})

test_that("skipUnusedFitPrep: a predict-only job schedules nothing from this module", {
  expect_identical(skipInit(list(doFit = FALSE)), character())
  expect_identical(skipInit(list(doFit = NA)), character())      # the ledger has this polygon: no fit
})

test_that("skipUnusedFitPrep: a job that fits, validates or holds out a fold still prepares", {
  expect_identical(skipInit(list(doFit = TRUE)), c("spreadFitPrepare", "run", "postFitDiagnostics"))
  expect_identical(skipInit(list(doFit = NA), list(studyAreaWithSpreadParams = toyLedger("9.8"))),
                   c("spreadFitPrepare", "run", "postFitDiagnostics"))          # no row for this polygon: fit
  expect_identical(skipInit(list(doFit = FALSE, mode = c("fit", "validate"))), c("spreadFitPrepare", "crossValidate"))
  expect_identical(skipInit(list(doFit = FALSE, heldOutFold = 1L)), c("spreadFitPrepare", "crossValidate"))
  expect_identical(skipInit(list(doFit = FALSE, mode = "profile")), c("spreadFitPrepare", "profileFromLedger"))   # a profile reads the covariates
})
