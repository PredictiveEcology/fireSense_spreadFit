## The final DEoptim population is written (DEoptimPopulation_<runName>.csv), and mode "terms" evaluates a pool of
## parameter vectors on the fit's data and writes the objective's separate terms (spreadFitTerms_<runName>.csv).
## See R/spreadFitPopulation.R.

modOut <- function(sim) file.path(SpaDES.core::outputPath(sim), toyModule)
termCols <- c("objective", "adTest", "pitAD", "yearArea", "areaDist", "crpsFire", "crpsTotalArea", "crpsTop10",
              "crpsTopShare", "nFires")

test_that("a fit writes its final population, one row per member: the parameters and the objective value", {
  sim <- toySim()
  mockFitAndLedger(sim)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  out <- list(sim = sim)
  f <- file.path(modOut(out$sim), "DEoptimPopulation_toyRun.csv")
  expect_true(file.exists(f))
  pop <- data.table::fread(f)
  parNames <- names(SpaDES.core::P(out$sim, module = toyModule)$lower)
  expect_identical(names(pop), c(parNames, "objFunVal"))
  expect_identical(nrow(pop), 7L)
  expect_equal(pop[[parNames[1]]], 1:7)                  # toyDE()'s final population
  expect_equal(pop$objFunVal, 7:1)
})

## The real objective simulates fires and belongs to fireSenseUtils (tested there, with scorePopulationTerms()), as the
## profile's is in test-profileMode.R: here scorePopulationTerms() is a recorder that returns terms from the parameters
## and the seed alone, as the real one does.
termsRun <- function(params = list(), reps = 1L, ledgerFn = profileLedger) {
  rec <- new.env()
  pars <- list(mode = "terms", termsReps = reps, cores = NA, doFit = TRUE)
  pars[names(params)] <- params
  sim <- toySim(pars)
  sim <- runEvents(runEvents(sim, "init"), "spreadFitPrepare")
  parNames <- names(SpaDES.core::P(sim, module = toyModule)$lower)
  best <- stats::setNames(seq_along(parNames) / 10, parNames)
  sim$studyAreaWithSpreadParams <- ledgerFn(parNames, best)
  mockInModule(sim, runDEoptim = function(...) stop("runDEoptim was called"),
               CacheGeo = function(...) stop("the ledger was written"))
  testthat::local_mocked_bindings(.package = "fireSenseUtils",
    scorePopulationTerms = function(pop, fn, reps, cl, seed, fnArgs) {
      rec$args <- list(pop = pop, reps = reps, cl = cl, seed = seed, fnArgs = fnArgs, fn = fn)
      jobs <- data.table::CJ(member = seq_len(nrow(pop)), rep = seq_len(reps))
      terms <- data.table::as.data.table(matrix(0, nrow(jobs), length(termCols), dimnames = list(NULL, termCols)))
      terms[, objective := vapply(seq_len(nrow(jobs)), function(i)
        sum(unlist(pop[jobs$member[i], ])) + seed + jobs$rep[i], numeric(1))]
      cbind(jobs, terms)
    })
  list(sim = runEvents(sim, "termsFromPool"), parNames = parNames, best = best, rec = rec)
}

test_that("mode 'terms' schedules only spreadFitPrepare and termsFromPool, whatever doFit says", {
  for (df in list(TRUE, NA)) {
    sim <- runEvents(toySim(list(mode = "terms", doFit = df)), "init")
    expect_identical(queued(sim), c("spreadFitPrepare", "termsFromPool"))
  }
})

test_that("mode 'terms' scores the ledger's sets, one row per member with every term, in a file", {
  out <- termsRun()
  tr <- out$sim$spreadFitTerms
  expect_identical(nrow(tr), 2L)                              # profileLedger() holds two sets
  expect_identical(names(tr), c("member", "rep", out$parNames, termCols))
  expect_equal(unlist(tr[1, out$parNames, with = FALSE]), out$best)
  expect_null(out$rec$args$cl)                                # cores = NA: evaluated here, data in the args
  expect_true(all(emulatorDataArgs %in% names(out$rec$args$fnArgs)))
  expect_identical(out$rec$args$fn, fireSenseUtils::.objfunSpreadFit)
  f <- file.path(modOut(out$sim), "spreadFitTerms_toyRun.csv")
  expect_true(file.exists(f))
  expect_equal(data.table::fread(f)$objective, tr$objective)
})

test_that("mode 'terms' gives scorePopulationTerms() the seed and reps, so a fixed seed repeats", {
  a <- termsRun(list(termsSeed = 3L), reps = 2L)
  b <- termsRun(list(termsSeed = 3L), reps = 2L)
  d <- termsRun(list(termsSeed = 4L), reps = 2L)
  expect_identical(c(a$rec$args$seed, a$rec$args$reps), c(3L, 2L))
  expect_identical(nrow(a$sim$spreadFitTerms), 4L)            # member x rep
  expect_equal(a$sim$spreadFitTerms, b$sim$spreadFitTerms)
  expect_false(isTRUE(all.equal(a$sim$spreadFitTerms$objective, d$sim$spreadFitTerms$objective)))
})

test_that("mode 'terms' stops when the ledger row is not this run's fit", {
  expect_error(termsRun(ledgerFn = function(p, b) { l <- profileLedger(p, b); l$polygonID <- "9.8"; l }),
               "ledger has no row for 9.9")
  expect_error(termsRun(ledgerFn = function(p, b) profileLedger(sub("class1", "classX", p), b)),
               "fitted with another formula")
})

test_that("mode 'terms' scores the sets of a csv: the population file, extra columns ignored", {
  parNames <- names(SpaDES.core::P(runEvents(runEvents(toySim(), "init"), "spreadFitPrepare"), module = toyModule)$lower)
  f <- withr::local_tempfile(fileext = ".csv")
  pop <- data.table::as.data.table(rbind(seq_along(parNames) / 10, seq_along(parNames) / 5, seq_along(parNames) / 2))
  data.table::setnames(pop, parNames)
  pop[, objFunVal := c(9, 8, 7)]
  data.table::fwrite(pop, f)
  out <- termsRun(list(termsPool = f))
  expect_identical(nrow(out$sim$spreadFitTerms), 3L)
  expect_equal(out$sim$spreadFitTerms[[parNames[2]]], pop[[parNames[2]]])
  ## a pool missing a parameter is refused
  data.table::fwrite(pop[, -1], f)
  expect_error(termsRun(list(termsPool = f)), "does not hold")
})
