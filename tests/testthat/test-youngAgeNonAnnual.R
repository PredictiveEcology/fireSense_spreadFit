## youngAge must stay mutually exclusive with every other non-annual covariate even when youngAge
## itself is a non-annual column (true for ELFs where youngAge comes from cohort fuel classes, not
## an annual layer). test-spreadFitPrep.R's fixture always has youngAge annual (helper-toyInputs.R),
## which is why fireSense_spreadFit.R:519-526 appending every non-annual column name to youngAge's
## own mutuallyExclusiveCols entry -- including "youngAge" itself -- was never caught: with youngAge
## annual, it is never one of the non-annual names appended.

P1 <- function(sim) SpaDES.core::params(sim)[[moduleName]]

## youngAge moved from the annual covariates into the non-annual table
youngAgeNonAnnualObjects <- function() {
  objs <- toyObjects()
  objs$fireSense_annualSpreadFitCovariates <- lapply(
    objs$fireSense_annualSpreadFitCovariates,
    function(dt) { dt[, youngAge := NULL]; dt }
  )
  ## same young/not-young pixels as the original annual youngAge (see toyObjects()):
  ## pixelID 2,3,4,55,56,57,88,89,90 -> 0,1,0,0,0,1,0,0,1
  objs$fireSense_nonAnnualSpreadFitCovariates[[1]][, youngAge := c(0, 1, 0, 0, 0, 1, 0, 0, 1)]
  objs
}

test_that("the default mutuallyExclusiveCols does not include youngAge itself", {
  objs <- youngAgeNonAnnualObjects()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE), objs)
  sim <- runEvents(runEvents(sim, "init"), "spreadFitPrepare")
  mec <- P1(sim)$mutuallyExclusiveCols
  expect_false("youngAge" %in% mec$youngAge)
  ## the real bug: "youngAge" was appended because it is a non-annual column name
  expect_true(all(c("class1", "class2", "nf") %in% mec$youngAge))
})

test_that("after exclusivity, youngAge stays 1 on young pixels (objective's covariate path)", {
  objs <- youngAgeNonAnnualObjects()
  sim <- toySim(list(stopIfNoPreRunFit = FALSE), objs)
  rec <- mockFitAndLedger(sim)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  a <- rec$deArgs

  expect_false("youngAge" %in% a$mutuallyExclusive$youngAge)

  ## mirrors fireSenseUtils::spreadProbFromIntegerCovs()'s own merge of one year's annual
  ## covariates onto the (shared) non-annual table, before it applies makeMutuallyExclusive()
  nonAnnual <- data.table::as.data.table(a$nonAnnualDTx1000$year2001_year2002)
  annual <- data.table::as.data.table(a$annualDTx1000$year2001)
  merged <- nonAnnual[annual, on = "pixelID"]

  out <- fireSenseUtils::makeMutuallyExclusive(dt = data.table::copy(merged),
                                               mutuallyExclusiveCols = a$mutuallyExclusive)
  young <- out$youngAge != 0
  expect_true(any(young))
  expect_true(all(out$youngAge[young] != 0))     # youngAge itself is never zeroed
  expect_true(all(unlist(out[young, c("class1", "class2", "nf"), with = FALSE]) == 0))
})
