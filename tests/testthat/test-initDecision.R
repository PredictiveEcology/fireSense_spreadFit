## What `init` decides: fit, skip, or stop. Toy inputs; no event after `init` is run here.

fitEvents <- c("spreadFitPrepare", "spreadFitPrepare", "estimateThreshold", "run", "postFitDiagnostics")
noFitMsg <- "There is no pre-run SpreadFit \\(sim\\$studyAreaWithSpreadParams\\), but parameter `stopIfNoPreRunFit` is `TRUE`"

initOnly <- function(params = list(), objects = list()) runEvents(toySim(params, objects), "init")

test_that("no ledger row for this polygon: stopIfNoPreRunFit decides between stopping and fitting", {
  expect_error(initOnly(list(stopIfNoPreRunFit = TRUE)), noFitMsg)
  expect_identical(queued(initOnly(list(stopIfNoPreRunFit = FALSE))), fitEvents)
})

test_that("stopIfNoPreRunFit defaults to TRUE: a bare module stops rather than start a fit", {
  expect_error(initOnly(), noFitMsg)
})

test_that("a ledger holding only neighbours' rows is not a fit for this polygon", {
  objs <- list(studyAreaWithSpreadParams = toyLedger(c("9.8", "19.9")))   # this polygon is "9.9"
  expect_error(initOnly(list(stopIfNoPreRunFit = TRUE), objs), noFitMsg)
  expect_identical(queued(initOnly(list(stopIfNoPreRunFit = FALSE), objs)), fitEvents)
})

test_that("a ledger row for this polygon: nothing but the preparation is scheduled", {
  objs <- list(studyAreaWithSpreadParams = toyLedger(c("9.8", "9.9")))
  expect_identical(queued(initOnly(list(stopIfNoPreRunFit = TRUE), objs)), "spreadFitPrepare")
  expect_identical(queued(initOnly(list(stopIfNoPreRunFit = FALSE), objs)), "spreadFitPrepare")
})

test_that("the ledger is keyed on .ELFind, not on .runName", {
  ## .runName matches a row, .ELFind does not: no fit is found
  objs <- list(studyAreaWithSpreadParams = toyLedger("toyRun"))
  expect_identical(queued(initOnly(list(stopIfNoPreRunFit = FALSE), objs)), fitEvents)
  ## without .ELFind, .inputObjects copies .runName into it, and the row is found
  objs[".ELFind"] <- list(NULL)                    # toySim() drops an object set to NULL
  expect_identical(queued(initOnly(list(stopIfNoPreRunFit = FALSE), objs)), "spreadFitPrepare")
})

test_that("decision table: refitExisting x stopIfNoPreRunFit, when the ledger HAS this polygon", {
  ## `refitExisting = TRUE` is an explicit request to fit, so it OVERRIDES `stopIfNoPreRunFit`:
  ## a refit does NOT need `stopIfNoPreRunFit = FALSE` as well.
  objs <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  cell <- function(refit, stop) initOnly(list(refitExisting = refit, stopIfNoPreRunFit = stop), objs)
  expect_identical(queued(cell(refit = FALSE, stop = FALSE)), "spreadFitPrepare")   # skip
  expect_identical(queued(cell(refit = FALSE, stop = TRUE)), "spreadFitPrepare")    # skip
  expect_identical(queued(cell(refit = TRUE, stop = FALSE)), fitEvents)             # refit
  expect_identical(queued(cell(refit = TRUE, stop = TRUE)), fitEvents)              # refit (override)
})

test_that("decision table: refitExisting when the ledger lacks this polygon", {
  cell <- function(refit, stop) initOnly(list(refitExisting = refit, stopIfNoPreRunFit = stop))
  expect_identical(queued(cell(refit = FALSE, stop = FALSE)), fitEvents)
  expect_identical(queued(cell(refit = TRUE, stop = FALSE)), fitEvents)
  expect_error(cell(refit = FALSE, stop = TRUE), noFitMsg)
  expect_identical(queued(cell(refit = TRUE, stop = TRUE)), fitEvents)              # refit (override)
})

test_that("refitExisting overrides the DEFAULT stopIfNoPreRunFit, without setting it", {
  ## the user-facing promise: setting only `refitExisting = TRUE` is enough to get a refit
  objs <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  expect_identical(queued(initOnly(list(refitExisting = TRUE), objs)), fitEvents)
  expect_identical(queued(initOnly(list(), objs)), "spreadFitPrepare")  # default still skips
  expect_identical(queued(initOnly(list(refitExisting = TRUE))), fitEvents)  # and with no ledger row
})

test_that("stopIfNoPreRunFit still stops when refitExisting is not exactly TRUE", {
  expect_error(initOnly(list(refitExisting = FALSE, stopIfNoPreRunFit = TRUE)), noFitMsg)
  expect_error(initOnly(list(refitExisting = NA, stopIfNoPreRunFit = TRUE)), noFitMsg)
})

test_that("refitExisting must be exactly TRUE: NA or a string does not trigger a refit", {
  objs <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  expect_identical(queued(initOnly(list(refitExisting = NA, stopIfNoPreRunFit = FALSE), objs)),
                   "spreadFitPrepare")
})

test_that("mode picks the events that follow the threshold", {
  prep <- c("spreadFitPrepare", "spreadFitPrepare", "estimateThreshold")
  p <- function(mode) list(stopIfNoPreRunFit = FALSE, mode = mode)
  expect_identical(queued(initOnly(p("debug"))), c(prep, "debug"))
  expect_identical(queued(initOnly(p(c("fit", "visualize")))),
                   c(prep, "run", "postFitDiagnostics", "debug", "plot"))
  ## 'debug' wins over 'visualize': no DEoptim run, so nothing to visualize
  expect_identical(queued(initOnly(p(c("debug", "visualize")))), c(prep, "debug"))
})

test_that("mode 'validate' adds crossValidate, after the fit or, with a ledger row, on its own", {
  p <- list(stopIfNoPreRunFit = FALSE, mode = c("fit", "validate"))
  expect_identical(queued(initOnly(p)), c(fitEvents, "crossValidate"))
  objs <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  ## the ledger's fit is not refitted, but validation needs the threshold
  expect_identical(queued(initOnly(p, objs)), c("spreadFitPrepare", "estimateThreshold", "crossValidate"))
})

test_that("supplying parsKnown switches the module to debug mode", {
  sim <- initOnly(list(stopIfNoPreRunFit = FALSE, mode = "fit"),
                  list(parsKnown = c(0.26, 1, 2, 1, -1, 1, 1, 1)))
  expect_identical(SpaDES.core::params(sim)[[moduleName]]$mode, c("fit", "debug"))
  expect_identical(tail(queued(sim), 1), "debug")
})

test_that("events are scheduled at .runInitialTime", {
  sim <- toySim(list(stopIfNoPreRunFit = FALSE, .runInitialTime = 1))
  sim <- runEvents(sim, "init")
  ev <- SpaDES.core::events(sim)
  expect_identical(unique(ev$eventTime[ev$moduleName == moduleName]), 1)
})

test_that("heldOutFold NA (the default) leaves scheduling unchanged", {
  expect_identical(queued(initOnly(list(stopIfNoPreRunFit = FALSE))),
                   queued(initOnly(list(stopIfNoPreRunFit = FALSE, heldOutFold = NA))))
})

test_that("heldOutFold 1 or 2 schedules ONLY prepare, threshold and crossValidate", {
  ## no full fit ("run"), no postFitDiagnostics -- and it bypasses stopIfNoPreRunFit entirely,
  ## which would otherwise stop a bare module (see the first test in this file)
  expect_identical(queued(initOnly(list(heldOutFold = 1L))),
                   c("spreadFitPrepare", "estimateThreshold", "crossValidate"))
  expect_identical(queued(initOnly(list(heldOutFold = 2L))),
                   c("spreadFitPrepare", "estimateThreshold", "crossValidate"))
  ## still true with a ledger row present and refitExisting set: heldOutFold does not consult either
  objs <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  expect_identical(queued(initOnly(list(heldOutFold = 1L, refitExisting = TRUE), objs)),
                   c("spreadFitPrepare", "estimateThreshold", "crossValidate"))
})

test_that("heldOutFold rejects values other than NA, 1 or 2", {
  expect_error(initOnly(list(heldOutFold = 3L)), "heldOutFold.*NA, 1L or 2L")
  expect_error(initOnly(list(heldOutFold = "a")))
})
