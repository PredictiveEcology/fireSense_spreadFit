## mode "profile": profile this ELF's ledger fit at given values (fireSenseUtils::profileCoefficients(at = )), no fit.
## See R/profileLedgerFit.R.

## A ledger row for the toy polygon with the toy formula's parameters, best set first
profileLedger <- function(parNames, best) {
  l <- toyLedger("9.9")
  sets <- data.table::as.data.table(rbind(best, best + 1))
  data.table::setnames(sets, parNames)
  l$params <- list(sets)
  l$objFunVal <- list(c(101, 102))
  l
}

profileRun <- function(params = list(), ledgerFn = NULL) {
  rec <- new.env()
  pars <- list(mode = "profile", profileAt = c(2, -5, 0), profileCoefs = "fuel", profileReps = 3L, cores = NA, doFit = TRUE)
  pars[names(params)] <- params
  sim <- toySim(pars)
  sim <- runEvents(runEvents(sim, "init"), "spreadFitPrepare") # the bounds, hence the parameter names
  parNames <- names(SpaDES.core::P(sim, module = toyModule)$lower)
  best <- stats::setNames(seq_along(parNames) / 10, parNames)
  sim$studyAreaWithSpreadParams <- if (is.null(ledgerFn)) profileLedger(parNames, best) else ledgerFn(parNames, best)
  mockInModule(sim, runDEoptim = function(...) stop("runDEoptim was called"),
               CacheGeo = function(...) stop("the ledger was written"))
  testthat::local_mocked_bindings(.package = "fireSenseUtils",
    profileCoefficients = function(best, pop, fn, reps, cl, fnArgs, coefs, at, ...) {
      rec$args <- list(best = best, pop = pop, reps = reps, cl = cl, fnArgs = fnArgs, coefs = coefs, at = at)
      data.table::data.table(coef = c("", coefs[1]), at = c(NA, at[[1]][1]), isZero = FALSE, mean = c(100, 110),
                             delta = c(0, 10), deltaSE = c(0, 1))
    })
  out <- runEvents(sim, "profileFromLedger")
  list(sim = out, rec = rec, best = best)
}

test_that("mode 'profile' schedules only spreadFitPrepare and profileFromLedger, whatever doFit says", {
  for (df in list(TRUE, NA)) {
    sim <- runEvents(toySim(list(mode = "profile", doFit = df)), "init")
    expect_identical(queued(sim), c("spreadFitPrepare", "profileFromLedger"))
  }
})

test_that("the profile is of the ledger's best set, at the given values, for the fuel coefficients only", {
  out <- profileRun()
  a <- out$rec$args
  expect_identical(a$best, out$best)
  expect_null(a$pop)
  expect_identical(a$coefs, c("class1", "class2"))                 # "fuel": the fuel biomass columns
  expect_identical(a$at, list(class1 = c(-5, 0, 2), class2 = c(-5, 0, 2)))
  expect_identical(a$reps, 3L)
  expect_null(a$cl)                                                # cores = NA: evaluated here, data in the args
  expect_true(all(emulatorDataArgs %in% names(a$fnArgs)))
  ## the result, with the ledger's own value of the best set, in the sim and in a file
  expect_identical(out$sim$spreadFitProfile$ledgerObjFunVal, c(101, 101))
  f <- file.path(SpaDES.core::outputPath(out$sim), toyModule, "spreadFitProfile_toyRun.csv")
  expect_true(file.exists(f))
  expect_equal(data.table::fread(f)$delta, c(0, 10))
})

test_that("profileCoefs NA profiles every covariate coefficient", {
  out <- profileRun(list(profileCoefs = NA_character_))
  expect_identical(out$rec$args$coefs, c("CMDsm", "youngAge", "class1", "class2", "nf"))
})

test_that("the profile stops when the ledger row is not this run's fit", {
  ## no row for this ELF
  expect_error(profileRun(ledgerFn = function(p, b) { l <- profileLedger(p, b); l$polygonID <- "9.8"; l }),
               "ledger has no row for 9.9")
  ## other parameters (another formula)
  expect_error(profileRun(ledgerFn = function(p, b) profileLedger(sub("class1", "classX", p), b)),
               "fitted with another formula")
  ## no values to profile at
  expect_error(profileRun(list(profileAt = NA_real_)), "needs `profileAt`")
})
