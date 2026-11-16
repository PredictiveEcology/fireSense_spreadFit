## After a fit that simulated its members, the fuel x climate interaction test of the fit's own fires
## (fireSenseUtils::fuelClimateInteractionTest()) goes in the ledger row, beside the in-sample sizes, and
## to a csv in the module's output folder. A fold, or a fit without simulations, has neither.

skip_if_not("fuelClimateInteractionTest" %in% getNamespaceExports("fireSenseUtils"))

## 8 fire years of 2 fires; the planted size depends on climate x fuel when `interaction` is not 0
toyInteractionCovs <- function(interaction) {
  withr::local_seed(2)
  years <- 2001:2008
  fuel <- rnorm(length(years) * 2L)
  clim <- rnorm(length(years))
  fy <- rep(seq_along(years), each = 2L)
  obs <- round(exp(5 + 0.3 * clim[fy] + 0.3 * fuel + interaction * clim[fy] * fuel + rnorm(length(fy), 0, 0.1)))
  nm <- paste0("year", years)
  pix <- seq_along(fy)
  list(
    covs = list(
      annualDTx1000 = stats::setNames(lapply(clim, function(cl) data.frame(pixelID = pix, CMDsm = cl * 1000)), nm),
      nonAnnualDTx1000 = list(`2000` = data.frame(pixelID = pix, class1 = fuel * 1000)),
      fireBufferedListDT = stats::setNames(lapply(seq_along(years), function(i)
        data.frame(pixelID = which(fy == i), buffer = 1L, ids = 1:2)), nm)),
    sims = data.table::data.table(member = 1L, yr = nm[fy], rep = 1L, ids = rep(1:2, length(years)),
                                  initialLocus = pix, size = obs, sim = obs))
}

test_that("interactionTest() tests the fit's own fires with the spread formula and the parameters' draws and seed", {
  x <- toyInteractionCovs(0.6)
  sim <- toySim(list(interactionTestDraws = 49L, interactionTestSeed = 7L),
                list(fireSense_spreadFormula = "~ 0 + CMDsm + class1"))
  DE <- toyDE(8L); attr(DE, "fitSims") <- x$sims
  r <- interactionTest(sim, DE, x$covs)
  expect_identical(r$nBoot, 49L)
  expect_identical(r$seed, 7L)
  expect_identical(r$climate, "CMDsm")
  t <- r$terms[r$terms$fuel == "class1" & r$terms$at == "burn" & r$terms$response == "logSize", ]
  expect_gt(t$bIntMedian, 0)
  expect_identical(r, interactionTest(sim, DE, x$covs))        # deterministic
  expect_null(interactionTest(sim, toyDE(8L), x$covs))         # no simulations
  ## a failure of the test is a warning, not the end of a fit
  bad <- x$covs; bad$fireBufferedListDT <- NULL
  expect_warning(expect_null(interactionTest(sim, DE, bad)), "interaction test was skipped")
  expect_null(interactionTest(toySim(list(interactionTestDraws = 0L)), DE, x$covs))
})

test_that("the run event writes the ledger column and the csv; a fit without the test has no column", {
  col <- fireSenseUtils::spreadFitInteractionTestTxt
  res <- list(terms = data.table::data.table(fuel = "class1", at = "burn", response = "logSize", bIntMedian = 0.5),
              climate = "CMDsm", nFires = 8L, nYears = 8L, nBoot = 999L, seed = 1L)
  fitWith <- function(test) {
    rec <- new.env()
    sim <- toySim()
    mockFitAndLedger(sim, rec)
    mockInModule(sim, interactionTest = function(sim, DE, covs) test)
    list(sim = suppressMessages(SpaDES.core::spades(sim)), row = rec$geoArgs$studyAreaFireSense)
  }
  unlink(list.files(file.path(toyPaths()$outputPath, toyModule), "^fuelClimateInteraction", full.names = TRUE))
  with <- fitWith(res)
  sim <- with$sim; row <- with$row
  expect_identical(row[[col]][[1]], res)
  csv <- data.table::fread(file.path(SpaDES.core::outputPath(sim), toyModule, "fuelClimateInteraction_toyRun.csv"))
  expect_identical(names(csv), c("run", "climate", "fuel", "at", "response", "bIntMedian"))
  expect_identical(csv$run, "toyRun")
  ## the row of a fit with no test (NULL) is what it was, and the old and the new bind
  without <- fitWith(NULL)$row
  expect_false(col %in% names(without))
  both <- dplyr::bind_rows(without, row)
  expect_identical(nrow(both), 2L)
  expect_true(is.null(both[[col]][[1]]) || all(is.na(both[[col]][[1]])))
})

test_that("with no simulations (simulateMembers = 0) the run event writes no test", {
  unlink(list.files(file.path(toyPaths()$outputPath, toyModule), "^fuelClimateInteraction", full.names = TRUE))
  rec <- new.env()
  sim <- toySim()
  mockFitAndLedger(sim, rec)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  expect_false(fireSenseUtils::spreadFitInteractionTestTxt %in% names(rec$geoArgs$studyAreaFireSense))
  expect_length(list.files(file.path(SpaDES.core::outputPath(sim), toyModule), pattern = "^fuelClimateInteraction"), 0L)
})
