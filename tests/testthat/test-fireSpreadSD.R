## fireSpreadSD, the sd of a per-fire random effect on spread (fireSenseUtils#221): fitted only with
## wholeLandscapeSpread, just before yearSpreadSD. The ledger row holds it and the escape definition.

wholeP <- list(wholeLandscapeSpread = TRUE, spreadLink = "powerLaw", powerLawGamma = 1.9, penaliseRunaways = FALSE)
moduleParamDefault <- function(nm) {
  p <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)$parameters
  p$default[[match(nm, p$paramName)]]
}
P1 <- function(sim) SpaDES.core::params(sim)[[moduleName]]
prepared <- function(params = list()) {
  runEvents(runEvents(toySim(params), "init"), "spreadFitPrepare")
}

test_that("estimateSpreadParams puts fireSpreadSD just before yearSpreadSD, and leaves it out unless given", {
  annual <- toyObjects()$fireSense_annualSpreadFitCovariates
  form <- "~ 0 + CMDsm + nf"
  up <- estimateSpreadParams(form, annual, "upper", 9, fireSpreadSDBounds = c(0, 10), yearSpreadSDBounds = c(0, 5))
  lo <- estimateSpreadParams(form, annual, "lower", 9, fireSpreadSDBounds = c(0, 10), yearSpreadSDBounds = c(0, 5))
  expect_identical(tail(names(up), 2), c("fireSpreadSD", "yearSpreadSD"))
  expect_identical(unname(c(lo["fireSpreadSD"], up["fireSpreadSD"])), c(0, 10))
  expect_identical(unname(c(lo["yearSpreadSD"], up["yearSpreadSD"])), c(0, 5))
  expect_identical(tail(names(estimateSpreadParams(form, annual, "upper", 9, fireSpreadSDBounds = c(0, 10))), 1), "fireSpreadSD")
  expect_false("fireSpreadSD" %in% names(estimateSpreadParams(form, annual, "upper", 9, yearSpreadSDBounds = c(0, 5))))
})

test_that("the default bounds have fireSpreadSD (0, 10) before yearSpreadSD when the fit is whole, and not otherwise", {
  expect_identical(moduleParamDefault("fireSpreadSDBounds"), c(0, 10))
  whole <- P1(prepared(wholeP))
  expect_identical(tail(names(whole$upper), 2), c("fireSpreadSD", "yearSpreadSD"))
  expect_identical(tail(names(whole$lower), 2), c("fireSpreadSD", "yearSpreadSD"))
  expect_identical(unname(c(whole$lower["fireSpreadSD"], whole$upper["fireSpreadSD"])), c(0, 10))
  notWhole <- P1(prepared())
  expect_false("fireSpreadSD" %in% c(names(notWhole$upper), names(notWhole$lower)))
  expect_identical(tail(names(notWhole$upper), 1), "yearSpreadSD")
})

test_that("a fit that is not whole says it leaves fireSpreadSD out", {
  expect_message(sd <- randomEffectBounds(list(fireSpreadSDBounds = c(0, 10), yearSpreadSDBounds = c(0, 10),
                                               wholeLandscapeSpread = FALSE)), "needs wholeLandscapeSpread")
  expect_identical(names(Filter(Negate(is.null), sd)), "yearSpreadSD")
  expect_silent(sd <- randomEffectBounds(list(fireSpreadSDBounds = c(0, 10), yearSpreadSDBounds = c(0, 10),
                                              wholeLandscapeSpread = TRUE)))
  expect_identical(names(Filter(Negate(is.null), sd)), c("fireSpreadSD", "yearSpreadSD"))
})

test_that("fireSpreadSDBounds = NA turns it off; a supplied bound without it keeps the other bound aligned", {
  off <- P1(prepared(c(wholeP, list(fireSpreadSDBounds = NA))))
  expect_false("fireSpreadSD" %in% c(names(off$upper), names(off$lower)))
  up <- c(maxAsymptote = 0.3, CMDsm = 1, youngAge = 0, class1 = 2, class2 = 3, nf = 4)
  p <- P1(prepared(c(wholeP, list(upper = up))))
  expect_identical(names(p$lower), names(up))
})

test_that("the ledger row holds fireSpreadSD as a fitted value, and the escape definition", {
  rec <- seededRun(params = wholeP)
  expect_true("fireSpreadSD" %in% names(rec$lower))
  row <- rec$geoArgs$studyAreaFireSense
  expect_true("fireSpreadSD" %in% names(row$params[[1]]))
  def <- fireSenseUtils::escapeDefinition()
  expect_identical(row$escapeSizeHa, def$escapeSizeHa)
  expect_identical(row$fireCauses, def$fireCauses)
})

test_that("the ledger row records the module's escapeSizeHa and fireCauses", {
  rec <- seededRun(params = list(escapeSizeHa = 20, fireCauses = "natural"))
  row <- rec$geoArgs$studyAreaFireSense
  expect_identical(row$escapeSizeHa, 20)
  expect_identical(row$fireCauses, "natural")
  expect_false("fireSpreadSD" %in% names(row$params[[1]]))
})

test_that("escapeSizeHa defaults to fireSenseUtils' constant, not a literal", {
  expect_identical(moduleParamDefault("escapeSizeHa"), fireSenseUtils::fireSenseEscapeSizeHa)
  expect_identical(moduleParamDefault("fireCauses"), fireSenseUtils::fireSenseFireCauses)
})

test_that("an emulator start has one column per bound, fireSpreadSD included, and the objective is told it is fitted", {
  rec <- seededRun(params = c(wholeP, list(initialPopMethod = "emulator")), mocks = emulatorMock)
  expect_identical(colnames(rec$emuPop), names(rec$lower))
  expect_identical(ncol(rec$calls[[1]]$DEoptimControl$initialpop), length(rec$lower))
  expect_true("fireSpreadSD" %in% colnames(rec$emuPop))
  expect_true(rec$emu$objArgs$fitFireSpreadSD)
  expect_true(rec$emu$objArgs$fitYearSpreadSD)
  expect_true(rec$emu$objArgs$wholeLandscape)
  ## not whole: nothing new reaches the objective
  rec2 <- seededRun(params = list(initialPopMethod = "emulator"), mocks = emulatorMock)
  expect_null(rec2$emu$objArgs$fitFireSpreadSD)
})
