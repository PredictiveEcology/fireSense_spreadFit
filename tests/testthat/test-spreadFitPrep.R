## The `spreadFitPrepare` event (spreadFitPrep()) on toy inputs. `mod$covsX1000` is not reachable
## from outside an event, so it is read from the arguments the (mocked) runDEoptim() receives in
## test-runEvent.R; this file checks what the event leaves in the simList.

prepared <- function(params = list(), objects = list()) {
  sim <- toySim(c(list(stopIfNoPreRunFit = FALSE), params), objects)
  runEvents(runEvents(sim, "init"), "spreadFitPrepare")
}
P1 <- function(sim) SpaDES.core::params(sim)[[moduleName]]

test_that("default bounds are built from the formula and the annual covariates, with no hillSlope1 or inflectionPoint1", {
  ## hillSlope1 and inflectionPoint1 are fixed at 1, not fitted (see NEWS and
  ## fireSenseUtils::fixLogisticPars()).
  p <- P1(prepared())
  expect_identical(p$upper, c(maxAsymptote = 0.276,
                              CMDsm = 50, youngAge = 0, class1 = 100, class2 = 100, nf = 50, yearSpreadSD = 1))
  expect_identical(p$lower, c(maxAsymptote = 0.25,
                              CMDsm = 0, youngAge = -50, class1 = -100, class2 = -100, nf = -50, yearSpreadSD = 0))
  expect_false(any(c("hillSlope1", "inflectionPoint1") %in% c(names(p$upper), names(p$lower))))
})

test_that("upperAndLowerValFuel sets the fuel bounds, and only those", {
  ## fuel is biomass / 1e4, not [0, 1]: bounds should be wide enough not to bind, not narrow
  p <- P1(prepared(list(upperAndLowerValFuel = 25)))
  expect_identical(unname(p$upper[c("class1", "class2", "nf", "CMDsm")]), c(25, 25, 50, 50))
  expect_identical(unname(p$lower[c("class1", "class2", "nf", "youngAge")]), c(-25, -25, -50, -50))
})

test_that("upperAndLowerVal sets the size of the default bounds", {
  p <- P1(prepared(list(upperAndLowerVal = 4)))
  expect_identical(unname(p$upper[c("CMDsm", "nf")]), c(4, 4))
  expect_identical(unname(p$upper["class1"]), 100)   # fuel has its own parameter
  expect_identical(unname(p$lower[c("youngAge", "nf")]), c(-4, -4))
})

test_that("supplied bounds are kept; only the missing one is filled", {
  up <- c(maxAsymptote = 0.3,
          CMDsm = 1, youngAge = 0, class1 = 2, class2 = 3, nf = 4)
  p <- P1(prepared(list(upper = up)))
  expect_identical(p$upper, up)
  expect_identical(unname(p$lower[c("class1", "nf")]), c(-100, -50))
  ## the supplied upper has no yearSpreadSD, so the filled-in lower has none either
  expect_identical(names(p$lower), names(up))
})

test_that("a supplied 'upper' or 'lower' naming hillSlope1 or inflectionPoint1 is refused", {
  up <- c(maxAsymptote = 0.3, hillSlope1 = 3, inflectionPoint1 = 5,
          CMDsm = 1, youngAge = 0, class1 = 2, class2 = 3, nf = 4)
  expect_error(prepared(list(upper = up)),
               "'hillSlope1', 'inflectionPoint1' found in the supplied 'upper'/'lower'")
  lo <- c(maxAsymptote = 0.2, hillSlope1 = 0.5, inflectionPoint1 = 0.1,
          CMDsm = 0, youngAge = -1, class1 = -2, class2 = -3, nf = -4)
  expect_error(prepared(list(lower = lo)),
               "'hillSlope1', 'inflectionPoint1' found in the supplied 'upper'/'lower'")
  ## each one on its own is refused too, and named
  expect_error(prepared(list(upper = up[names(up) != "hillSlope1"])),
               "'inflectionPoint1' found in the supplied 'upper'/'lower'")
  expect_error(prepared(list(lower = lo[names(lo) != "inflectionPoint1"])),
               "'hillSlope1' found in the supplied 'upper'/'lower'")
})

test_that("yearSpreadSD is in both bounds, last, unless turned off", {
  p <- P1(prepared())
  expect_identical(names(p$upper)[length(p$upper)], "yearSpreadSD")
  expect_identical(unname(c(p$lower["yearSpreadSD"], p$upper["yearSpreadSD"])), c(0, 1))
  off <- P1(prepared(list(yearSpreadSDBounds = NA)))
  expect_false("yearSpreadSD" %in% c(names(off$upper), names(off$lower)))
})

test_that("bounds whose names differ in order are refused", {
  up <- c(CMDsm = 1, maxAsymptote = 0.3,
          youngAge = 0, class1 = 2, class2 = 3, nf = 4)
  expect_error(prepared(list(upper = up)),
               "please ensure 'upper' and 'lower' params are named with an identical order")
})

test_that("the default mutuallyExclusiveCols gains every non-annual covariate; a custom one is kept", {
  expect_identical(P1(prepared())$mutuallyExclusiveCols,
                   list(youngAge = c("class", "nonForest", "class1", "class2", "nf")))
  custom <- list(youngAge = "class1")
  expect_identical(P1(prepared(list(mutuallyExclusiveCols = custom)))$mutuallyExclusiveCols, custom)
})

test_that("covMinMax_spread: fixed range for fuel, own range for cover and annual covariates", {
  sim <- prepared()
  ## fuel is biomass / 1e4 whatever the data hold (class1 reaches 22601); nf 0-0.8; CMDsm is CMDsm / 100 by default
  ## (covFixedRange), not the 10-40 these data span
  expect_identical(as.list(sim$covMinMax_spread),
                   list(class1 = c(0, 1e4), class2 = c(0, 1e4), nf = c(0, 0.8),
                        CMDsm = c(0, 100), youngAge = c(0, 1)))
})

test_that("covFixedRange = list() goes back to the data's range for CMDsm", {
  expect_identical(prepared(list(covFixedRange = list()))$covMinMax_spread$CMDsm, c(10, 40))
})

test_that("rescaleAll = FALSE leaves covMinMax_spread unset", {
  expect_null(prepared(list(rescaleAll = FALSE))$covMinMax_spread)
})

test_that("lociList: one table per year, fire size in pixels, ignition cell from the point", {
  ll <- prepared()$lociList
  expect_identical(names(ll), c("year2001", "year2002"))
  ## 5 ha and 12 ha on 1-ha pixels; (250, 950) is cell 3 and (550, 450) is row 6, column 6 = cell 56
  expect_equal(as.data.frame(ll$year2001),
               data.frame(size = c(5L, 12L), date = "year2001", ids = 11:12, cells = c(3L, 56L)))
  ## (850, 150) is row 9, column 9 = cell 89
  expect_equal(as.data.frame(ll$year2002),
               data.frame(size = 3L, date = "year2002", ids = 21L, cells = 89L))
})

test_that("fire size is in pixels of the raster, not hectares", {
  ## same fires on 200 m pixels (4 ha each): 12 ha = 3 pixels, 5 ha = 1.25 -> 1, 3 ha = 0.75 -> 1
  rtm <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, crs = toyCRS, vals = 1L)
  ll <- prepared(objects = list(rasterToMatch = rtm))$lociList
  expect_identical(ll$year2001$size, c(1L, 3L))
  expect_identical(ll$year2001$cells, c(2L, 13L))     # (250, 950) -> row 1, col 2; (550, 450) -> row 3, col 3
  expect_identical(ll$year2002$size, 1L)
})

test_that("a non-annual covariate that is all zeros is refused", {
  objs <- toyObjects()
  objs$fireSense_nonAnnualSpreadFitCovariates[[1]][, nf := 0]
  expect_error(prepared(objects = objs["fireSense_nonAnnualSpreadFitCovariates"]),
               "each non-annual spreadFit covariate cannot be all zeros")
})

test_that("a negative trace is refused", {
  expect_error(prepared(list(trace = -1)), "parameter 'trace' must be postive")
})

test_that("an NA in a covariate stops the fit, naming covMinMax_spread", {
  objs <- toyObjects()
  objs$fireSense_annualSpreadFitCovariates$year2001[2, CMDsm := NA]
  expect_error(prepared(objects = objs["fireSense_annualSpreadFitCovariates"]),
               "covMinMax_spread contains NA values")
})
