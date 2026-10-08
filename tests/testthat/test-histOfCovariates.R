## histOfCovariates(): what goes into the two covariate histograms. The plots are not drawn; the data
## and mappings they hold are checked.

dt <- data.table::data.table
annual <- list(year2001 = dt(pixelID = 1:3, CMDsm = c(10, 20, 30)),
               year2002 = dt(pixelID = 1:2, CMDsm = c(15, 25)))
## non-annual covariates are log biomass: log(1), log(2), log(10), log(4)
nonAnnual <- list(year2001_year2002 = dt(pixelID = 1:2, class1 = log(c(1, 2)), class2 = log(c(10, 4))))

hists <- function(...) {
  grDevices::pdf(NULL)                 # the function calls par(), which opens a device
  on.exit(grDevices::dev.off())
  histOfCovariates(...)
}

test_that("the annual plot holds every annual row, labelled by year, and drops pixelID", {
  out <- hists(annualList = annual, nonAnnualList = nonAnnual)
  expect_identical(names(out), c("annual", "nonAnnual"))
  d <- out$annual$data
  expect_identical(names(d), c("year", "Covariate", "value"))
  expect_identical(d$year, rep(c("year2001", "year2002"), c(3, 2)))
  expect_identical(as.character(d$Covariate), rep("CMDsm", 5))
  expect_identical(d$value, c(10, 20, 30, 15, 25))
})

## real ELFs use CMD, CMD_sp, cumMDC-derived columns, etc., not CMDsm; a plot hard-coded to CMDsm
## built without error but failed as soon as it was drawn (as `Plots()` does), with
## "object 'CMDsm' not found", stopping spreadFitPrepare() for every one of them
annualCMD <- list(year2001 = dt(pixelID = 1:3, CMD = c(10, 20, 30)),
                   year2002 = dt(pixelID = 1:2, CMD = c(15, 25)))

test_that("the annual plot renders for a climate covariate other than CMDsm", {
  out <- hists(annualList = annualCMD, nonAnnualList = nonAnnual)
  expect_no_error(ggplot2::ggplot_build(out$annual))
})

test_that("the annual plot works whatever the climate covariate is named, not just CMDsm", {
  out <- hists(annualList = annualCMD, nonAnnualList = nonAnnual)
  d <- out$annual$data
  expect_identical(names(d), c("year", "Covariate", "value"))
  expect_identical(as.character(d$Covariate), rep("CMD", 5))
  expect_identical(d$value, c(10, 20, 30, 15, 25))
})

test_that("the annual plot facets by covariate when there is more than one annual column", {
  annual2 <- list(year2001 = dt(pixelID = 1:2, CMD = c(10, 20), youngAge = c(1, 2)),
                  year2002 = dt(pixelID = 1:2, CMD = c(15, 25), youngAge = c(3, 4)))
  out <- hists(annualList = annual2, nonAnnualList = nonAnnual)
  d <- out$annual$data
  expect_identical(names(d), c("year", "Covariate", "value"))
  expect_setequal(unique(as.character(d$Covariate)), c("CMD", "youngAge"))
  expect_identical(d$value[as.character(d$Covariate) == "CMD"], c(10, 20, 15, 25))
  expect_identical(d$value[as.character(d$Covariate) == "youngAge"], c(1, 2, 3, 4))
})

test_that("the non-annual plot is long format, back-transformed from log, one row per pixel x fuel", {
  out <- hists(annualList = annual, nonAnnualList = nonAnnual)
  d <- out$nonAnnual$data
  expect_identical(names(d), c("year", "Fuel", "LogBiomass"))
  expect_identical(as.character(d$Fuel), rep(c("class1", "class2"), each = 2))
  ## exp(log(x)) = x; the minimum (log(1) = 0 -> 1) is set to exactly 1
  expect_equal(d$LogBiomass, c(1, 2, 10, 4))
  expect_identical(d$LogBiomass[1], 1)
})

test_that("the smallest value is forced to 1 whatever it back-transforms to", {
  ## minimum log value -2 -> exp(-2) = 0.135, replaced by 1 so the log10 axis has no sub-1 bar
  na2 <- list(a = dt(pixelID = 1:3, class1 = c(-2, 0, log(5))))
  d <- hists(annualList = annual, nonAnnualList = na2)$nonAnnual$data
  expect_equal(d$LogBiomass, c(1, 1, 5))
})

test_that("the inputs are not changed by reference", {
  a <- data.table::copy(annual); n <- data.table::copy(nonAnnual)
  hists(annualList = a, nonAnnualList = n)
  expect_identical(names(a$year2001), c("pixelID", "CMDsm"))
  expect_identical(n$year2001_year2002$class1, log(c(1, 2)))
})
