## resimulateCutFires: the whole-ELF covariate tables reach the objective on every path that calls it, and
## FALSE leaves the objective's arguments and the fit's cache key as they were.
## Helpers: toyWholeObjects(), fitArgs(), omitArgsCO() (helper-*.R).

wholeParams <- list(resimulateCutFires = TRUE, penaliseRunaways = FALSE, cutMaxTries = 20L)
newArgs <- c(wholeCovNames, "cutMaxTries")

test_that("FALSE passes the objective nothing new, and the fit's cache key is what it was", {
  expect_identical(wholeDeArgs(list(), list(resimulateCutFires = FALSE, cutMaxTries = 50L)), list())
  fit <- fitArgs(list(resimulateCutFires = FALSE), toyWholeObjects())
  expect_false(any(newArgs %in% names(fit)))
  covs <- list(annualDTx1000 = 1)
  expect_identical(c(list("run", objectiveBodies()), wholeCacheExtra(covs, wholeDeArgs(covs, list()))),
                   list("run", objectiveBodies()))
  omit <- omitArgsCO(list(resimulateCutFires = FALSE))
  expect_false(any(newArgs %in% unlist(omit$omit)))
})

test_that("TRUE gives runDEoptim() the whole tables and cutMaxTries, prepared as the buffer tables", {
  fit <- fitArgs(wholeParams, toyWholeObjects())
  expect_named(fit$annualDTx1000Whole, c("year2001", "year2002"))
  expect_identical(fit$cutMaxTries, 20L)
  ## same integers x 1000 as the buffer table on the buffer pixels, plus the pixel outside
  buf <- fit$annualDTx1000$year2001
  who <- fit$annualDTx1000Whole$year2001
  expect_identical(who$pixelID, c(1L, buf$pixelID))
  expect_identical(who$CMDsm[-1], buf$CMDsm)
  expect_identical(fit$covCentre, fitArgs(list(), list())$covCentre)
  ## fuel on the linear scale, as in the buffer table
  nb <- fit$nonAnnualDTx1000[[1]]; nw <- fit$nonAnnualDTx1000Whole[[1]]
  expect_identical(nw$class1[nw$pixelID %in% nb$pixelID], nb$class1)
})

test_that("TRUE keeps the big tables out of the fit's Cache key, and puts their digest in", {
  o <- toyWholeObjects() ## toyObjectsCO() names its data year "2000"
  o$fireSense_nonAnnualSpreadFitCovariatesWhole <- stats::setNames(o$fireSense_nonAnnualSpreadFitCovariatesWhole, "2000")
  omit <- omitArgsCO(wholeParams, o)
  expect_true(all(wholeCovNames %in% unlist(omit$omit)))
  covs <- list(wholeDigest = "abc")
  expect_identical(wholeCacheExtra(covs, wholeDeArgs(c(covs, list(annualDTx1000Whole = 1, nonAnnualDTx1000Whole = 2)),
                                                    list(resimulateCutFires = TRUE, cutMaxTries = 5L))),
                   list(wholeCovsDigest = "abc"))
})

test_that("TRUE with penaliseRunaways, or without the whole tables, stops", {
  expect_error(fitArgs(list(resimulateCutFires = TRUE, penaliseRunaways = TRUE), toyWholeObjects()),
               "cannot be combined with penaliseRunaways")
  expect_error(fitArgs(list(resimulateCutFires = TRUE, penaliseRunaways = FALSE)), "needs fireSense_annualSpreadFitCovariatesWhole")
})

test_that("whole tables with other years than the buffer tables stop", {
  o <- toyWholeObjects()
  o$fireSense_nonAnnualSpreadFitCovariatesWhole <- stats::setNames(o$fireSense_nonAnnualSpreadFitCovariatesWhole, "other")
  expect_error(fitArgs(wholeParams, o), "same data years")
  o <- toyWholeObjects()
  o$fireSense_annualSpreadFitCovariatesWhole <- o$fireSense_annualSpreadFitCovariatesWhole[1]
  expect_error(fitArgs(wholeParams, o), "lacks the fire years")
})

test_that("a cross-validation fold gets its own years of the whole table", {
  covs <- list(annualDTx1000 = list(y1 = 1, y2 = 2, y3 = 3, y4 = 4), fireBufferedListDT = list(y1 = 1, y2 = 2, y3 = 3, y4 = 4),
               historicalFires = list(y1 = 1, y2 = 2, y3 = 3, y4 = 4),
               annualDTx1000Whole = list(y1 = 1, y2 = 2, y3 = 3, y4 = 4), nonAnnualDTx1000Whole = list(a = 1), wholeDigest = "d")
  seen <- NULL
  local_mocked_bindings(fitSpread = function(sim, covs, ...) { seen <<- covs; stop("stop here") },
)
  sim <- toySim(wholeParams, toyWholeObjects())
  e <- tryCatch(fitAndScoreFold(sim, covs, fold = cvFolds(names(covs$historicalFires)), k = 1L), error = function(e) e)
  expect_named(seen$annualDTx1000Whole, names(seen$annualDTx1000))
  expect_identical(seen$nonAnnualDTx1000Whole, covs$nonAnnualDTx1000Whole)
})

test_that("the debug event and the emulator carry the whole tables", {
  expect_true(all(wholeCovNames %in% emulatorShipped(list(annualDTx1000Whole = 1, nonAnnualDTx1000Whole = 2))))
  expect_identical(emulatorShipped(list()), emulatorDataArgs)
  seen <- NULL
  sim <- toySim(c(list(mode = "debug"), wholeParams), toyWholeObjects())
  mockFitAndLedger(sim, new.env())
  mockInModule(sim, runSpreadWithoutDEoptim = function(...) { seen <<- list(...); NULL })
  suppressMessages(SpaDES.core::spades(sim))
  expect_true(all(newArgs %in% names(seen$objfunArgs)))
  expect_identical(seen$objfunArgs$cutMaxTries, 20L)
  ## FALSE: nothing extra
  sim <- toySim(c(list(mode = "debug"), list(resimulateCutFires = FALSE)), toyWholeObjects())
  mockFitAndLedger(sim, new.env())
  mockInModule(sim, runSpreadWithoutDEoptim = function(...) { seen <<- list(...); NULL })
  suppressMessages(SpaDES.core::spades(sim))
  expect_identical(seen$objfunArgs, list())
})

## what the emulator's Cache() is given, with the whole tables (or without)
emulatorCacheCall <- function(objects, params = wholeParams) {
  rec <- new.env()
  sim <- toySim(c(list(initialPopMethod = "emulator"), params), objects)
  mockFitAndLedger(sim, rec)
  mockInModule(sim, Cache = function(FUN, ..., .functionName = "", omitArgs = NULL, .cacheExtra = NULL) {
    if (startsWith(.functionName, "emulatorStartPopulation")) {
      rec$emu <- list(args = list(...), omitArgs = omitArgs, cacheExtra = .cacheExtra)
      lower <- rec$deLower <- list(...)$lower
      return(list(pop = matrix(rep(lower, each = 10L), 10L, dimnames = list(NULL, names(lower)))))
    }
    if (is.function(FUN)) FUN(...) else FUN
  })
  suppressMessages(SpaDES.core::spades(sim))
  rec$emu
}

test_that("the emulator's cache key follows the whole tables' digest, and never digests the tables", {
  a <- emulatorCacheCall(toyWholeObjects())
  expect_true("wholeTables" %in% a$omitArgs)
  expect_false(any(wholeCovNames %in% names(a$args$objArgs)))
  expect_true(all(wholeCovNames %in% names(a$args$wholeTables)))
  o <- toyWholeObjects()
  o$fireSense_annualSpreadFitCovariatesWhole$year2001$CMDsm[1] <- 99
  b <- emulatorCacheCall(o)
  expect_false(identical(a$cacheExtra, b$cacheExtra))
  expect_identical(a$cacheExtra[-length(a$cacheExtra)], b$cacheExtra[-length(b$cacheExtra)])
  ## FALSE: the key is what it was (no digest, no tables)
  off <- emulatorCacheCall(toyWholeObjects(), list(resimulateCutFires = FALSE))
  expect_length(off$cacheExtra, 2L)
  expect_length(off$args$wholeTables, 0L)
})
