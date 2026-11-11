## spreadLink = "climateOutside": the logistic and its maxAsymptote ceiling are fitted to the non-climate covariates
## and the climate term is added on the logit after it (fireSenseUtils::spreadProbFromCovs(climateRef =)). The
## reference climate is found once from the data the objective uses (fireSenseUtils::spreadClimateRef()) and handed
## to the fit, the debug objective call, the held-out simulation and the ledger row, which
## records it with the link. With the default "ceiling" nothing changes: no reference, no ledger columns, and the
## cache keys of the fit and the held-out simulation are what they were.

## the mean of the rescaled CMDsm over the 9 toy pixel-years, with the range the fit stores
cmdMean <- function(sim) {
  cmd <- c(10, 20, 30, 40, 25, 15, 12, 22, 32)
  rng <- sim$covMinMax_spread$CMDsm
  mean((cmd - rng[1]) / (rng[2] - rng[1]))
}

test_that("the 'ceiling' link leaves the fit, the ledger row and the cache keys as they were", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(md$parameters$default[[which(md$parameters$paramName == "spreadLink")]], "powerLaw") # the module's default
  expect_identical(fireSenseUtils::spreadLinkOptions[1], "ceiling")
  out <- fittedCO()
  expect_null(out$sim$climateRef_spread)
  expect_null(out$rec$deArgs$climateRef)                 # given to runDEoptim() as NULL, like covCentre
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_false(refCol %in% names(row))
  expect_false(linkCol %in% names(row))
  expect_identical(setdiff(names(row), ledgerKeyCols), fireSenseUtils::spreadFitAdditionalColNamesTxt)
})

test_that("with climateOutside the fit gets the climate reference, found once from the data the objective uses", {
  out <- fittedCO(list(spreadLink = "climateOutside"))
  ref <- out$sim$climateRef_spread
  expect_named(ref, "CMDsm")                       # the one climate covariate; youngAge is an indicator
  expect_equal(ref$CMDsm, cmdMean(out$sim))
  a <- out$rec$deArgs
  expect_identical(a$climateRef, ref)
  expect_equal(ref, fireSenseUtils::spreadClimateRef(a$annualDTx1000, a$nonAnnualDTx1000, a$formulaToFit,
                                                     a$covMinMax, a$mutuallyExclusive, covCentre = a$covCentre))
  ## the climate coefficient keeps its sign constraint
  expect_identical(unname(a$lower[["CMDsm"]]), 0)
})

test_that("the ledger row records the link and the reference; the parameter sets are unchanged in form", {
  out <- fittedCO(list(spreadLink = "climateOutside"))
  row <- out$rec$geoArgs$studyAreaFireSense
  expect_identical(row[[linkCol]], "climateOutside")
  expect_identical(row[[refCol]][[1]], out$sim$climateRef_spread)
  expect_identical(names(row$params[[1]])[1:4], c("maxAsymptote", "hillSlope1", "inflectionPoint1", "CMDsm"))
  expect_identical(out$sim$studyAreaWithSpreadParams, row)
})

test_that("with an intercept the reference is on the centred scale: 0", {
  out <- fittedCO(list(spreadLink = "climateOutside"), formula = sub("0 \\+", "1 +", formCO))
  expect_equal(out$sim$climateRef_spread$CMDsm, 0)
  expect_equal(out$sim$covCentre_spread$CMDsm, cmdMean(out$sim))
  expect_identical(out$rec$deArgs$climateRef, out$sim$climateRef_spread)
})

test_that("climateOutside without a climate covariate in the formula stops, saying so", {
  sim <- toySim(list(spreadLink = "climateOutside"), toyObjectsCO("~ 0 + youngAge + class1 + class2 + nf"))
  mockFitAndLedger(sim)
  expect_error(suppressMessages(SpaDES.core::spades(sim)), "needs a climate covariate")
})

test_that("the fit's cache key omits climateRef only when it is NULL", {
  off <- omitArgsCO()
  expect_length(off$omit, 1L)
  expect_true("climateRef" %in% off$omit[[1]])
  on <- omitArgsCO(list(spreadLink = "climateOutside"))
  expect_false("climateRef" %in% on$omit[[1]])
})

test_that("runSpreadWithoutDEoptim() hands climateRef to the objective, as it does covCentre", {
  src <- paste(deparse(body(runSpreadWithoutDEoptim)), collapse = "\n")
  expect_equal(lengths(regmatches(src, gregexpr("climateRef = climateRef", src))), 1L)
  expect_null(formals(runSpreadWithoutDEoptim)$climateRef)
})

test_that("the spread link is in the objective bodies the fit's cache key digests", {
  expect_true("spreadProbFromCovs" %in% names(objectiveBodies()))
})

test_that("the held-out simulation and the fold's ledger row use the fit's reference", {
  run <- function(params) {
    rec <- new.env(); rec$fnArgs <- list()
    sim <- toySim(c(list(heldOutFold = 1L, simulateMembers = 2L, doFit = FALSE), params), toyObjectsCO())
    mockFitAndLedger(sim, rec)
    mockInModule(sim, Cache = function(FUN, pop, fnArgs, ...) {
      if (!is.function(FUN)) return(FUN)
      rec$fnArgs <- fnArgs
      data.table::data.table(member = 1L, yr = "year2001", rep = 1L, ids = 1L, size = 4, sim = 8L)
    })
    list(sim = suppressMessages(SpaDES.core::spades(sim)), rec = rec)
  }
  on <- run(list(spreadLink = "climateOutside"))
  expect_identical(on$rec$fnArgs$climateRef, on$sim$climateRef_spread)
  expect_identical(on$sim$spreadFitHeldOut$fit[[refCol]][[1]], on$sim$climateRef_spread)
  expect_identical(on$sim$spreadFitHeldOut$fit[[linkCol]], "climateOutside")
  off <- run(list())
  expect_false("climateRef" %in% names(off$rec$fnArgs))
  expect_false(refCol %in% names(off$sim$spreadFitHeldOut$fit))
})
