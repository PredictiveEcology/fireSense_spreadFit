## What `init` decides: fit, skip, or stop. Toy inputs; no event after `init` is run here.

fitEvents <- c("spreadFitPrepare", "estimateThreshold", "run", "postFitDiagnostics")
noFitMsg <- "no usable pre-run fit in the ledger for polygon 9\\.9.*doFit"

## the messages `init` sends (runEvents() hides them)
initMessages <- function(params = list(), objects = list()) {
  msgs <- character()
  withCallingHandlers(
    try(SpaDES.core::spades(toySim(params, objects), events = stats::setNames(list("init"), moduleName)),
        silent = TRUE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  msgs
}

initOnly <- function(params = list(), objects = list()) runEvents(toySim(params, objects), "init")

test_that("doFit x ledger, when the ledger LACKS this polygon", {
  expect_identical(queued(initOnly(list(doFit = NA))), fitEvents)    # whatever is needed: fit
  expect_identical(queued(initOnly()), fitEvents)                    # NA is the default
  expect_identical(queued(initOnly(list(doFit = TRUE))), fitEvents)
  expect_error(initOnly(list(doFit = FALSE)), noFitMsg)             # never fit: stop, naming the polygon
})

test_that("a ledger holding only neighbours' rows is not a fit for this polygon", {
  objs <- list(studyAreaWithSpreadParams = toyLedger(c("9.8", "19.9")))   # this polygon is "9.9"
  expect_error(initOnly(list(doFit = FALSE), objs), noFitMsg)
  expect_identical(queued(initOnly(list(doFit = NA), objs)), fitEvents)
})

test_that("doFit x ledger, when the ledger HAS this polygon", {
  objs <- list(studyAreaWithSpreadParams = toyLedger(c("9.8", "9.9")))
  expect_identical(queued(initOnly(list(doFit = NA), objs)), "spreadFitPrepare")   # use the ledger fit
  expect_identical(queued(initOnly(list(), objs)), "spreadFitPrepare")
  expect_identical(queued(initOnly(list(doFit = FALSE), objs)), "spreadFitPrepare")
  expect_identical(queued(initOnly(list(doFit = TRUE), objs)), fitEvents)          # replace it
})

test_that("the ledger is keyed on .ELFind, not on .runName", {
  ## .runName matches a row, .ELFind does not: no fit is found
  objs <- list(studyAreaWithSpreadParams = toyLedger("toyRun"))
  expect_identical(queued(initOnly(list(), objs)), fitEvents)
  ## without .ELFind, .inputObjects copies .runName into it, and the row is found
  objs[".ELFind"] <- list(NULL)                    # toySim() drops an object set to NULL
  expect_identical(queued(initOnly(list(), objs)), "spreadFitPrepare")
})

test_that("refitExisting and stopIfNoPreRunFit are deprecated aliases of doFit, with one message", {
  objs <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  expect_length(grep("deprecated.*doFit", initMessages(list(refitExisting = TRUE), objs)), 1L)
  expect_identical(queued(initOnly(list(refitExisting = TRUE), objs)), fitEvents)                                        # refitExisting TRUE -> doFit TRUE
  expect_length(grep("deprecated.*doFit", initMessages(list(stopIfNoPreRunFit = TRUE))), 1L)
  expect_error(suppressMessages(initOnly(list(stopIfNoPreRunFit = TRUE))), noFitMsg)  # -> doFit FALSE
  ## stopIfNoPreRunFit = FALSE only said "allow a fit": that is what NA does
  expect_identical(queued(suppressMessages(initOnly(list(stopIfNoPreRunFit = FALSE)))), fitEvents)
  expect_identical(queued(suppressMessages(initOnly(list(stopIfNoPreRunFit = FALSE), objs))), "spreadFitPrepare")
  ## both set: refitExisting wins, as before
  expect_identical(queued(suppressMessages(initOnly(list(refitExisting = TRUE, stopIfNoPreRunFit = TRUE), objs))),
                   fitEvents)
  expect_identical(queued(suppressMessages(initOnly(list(refitExisting = FALSE, stopIfNoPreRunFit = TRUE), objs))),
                   "spreadFitPrepare")
  expect_error(suppressMessages(initOnly(list(refitExisting = FALSE, stopIfNoPreRunFit = TRUE))), noFitMsg)
  ## an explicit doFit is not overridden by an alias
  expect_identical(queued(suppressMessages(initOnly(list(doFit = FALSE, refitExisting = TRUE), objs))),
                   "spreadFitPrepare")
})

test_that("a run that sets none of the deprecated parameters prints no deprecation message", {
  expect_length(grep("deprecated", initMessages(list(doFit = NA))), 0L)
})

test_that("mode picks the events that follow the threshold", {
  prep <- c("spreadFitPrepare", "estimateThreshold")
  p <- function(mode) list(mode = mode)
  expect_identical(queued(initOnly(p("debug"))), c(prep, "debug"))
  expect_identical(queued(initOnly(p(c("fit", "visualize")))),
                   c(prep, "run", "postFitDiagnostics", "debug", "plot"))
  ## 'debug' wins over 'visualize': no DEoptim run, so nothing to visualize
  expect_identical(queued(initOnly(p(c("debug", "visualize")))), c(prep, "debug"))
})

test_that("mode 'validate' adds crossValidate, after the fit or, with a ledger row, on its own", {
  p <- list(mode = c("fit", "validate"))
  expect_identical(queued(initOnly(p)), c(fitEvents, "crossValidate"))
  objs <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  ## the ledger's fit is not refitted, but validation needs the threshold
  expect_identical(queued(initOnly(p, objs)), c("spreadFitPrepare", "estimateThreshold", "crossValidate"))
})

test_that("supplying parsKnown switches the module to debug mode", {
  sim <- initOnly(list(mode = "fit"),
                  list(parsKnown = c(0.26, 1, 2, 1, -1, 1, 1, 1)))
  expect_identical(SpaDES.core::params(sim)[[moduleName]]$mode, c("fit", "debug"))
  expect_identical(tail(queued(sim), 1), "debug")
})

test_that("events are scheduled at .runInitialTime", {
  sim <- toySim(list(.runInitialTime = 1))
  sim <- runEvents(sim, "init")
  ev <- SpaDES.core::events(sim)
  expect_identical(unique(ev$eventTime[ev$moduleName == moduleName]), 1)
})

test_that("heldOutFold NA (the default) leaves scheduling unchanged", {
  expect_identical(queued(initOnly()),
                   queued(initOnly(list(heldOutFold = NA))))
})

test_that("heldOutFold 1 or 2 schedules ONLY prepare, threshold and crossValidate", {
  ## no full fit ("run"), no postFitDiagnostics -- and it bypasses doFit entirely,
  ## so even doFit = FALSE without a ledger row does not stop it
  expect_identical(queued(initOnly(list(heldOutFold = 1L, doFit = FALSE))),
                   c("spreadFitPrepare", "estimateThreshold", "crossValidate"))
  expect_identical(queued(initOnly(list(heldOutFold = 2L))),
                   c("spreadFitPrepare", "estimateThreshold", "crossValidate"))
  ## still true with a ledger row present and doFit = TRUE: heldOutFold does not consult it
  objs <- list(studyAreaWithSpreadParams = toyLedger("9.9"))
  expect_identical(queued(initOnly(list(heldOutFold = 1L, doFit = TRUE), objs)),
                   c("spreadFitPrepare", "estimateThreshold", "crossValidate"))
})

test_that("heldOutFold rejects values other than NA, 1 or 2", {
  expect_error(initOnly(list(heldOutFold = 3L)), "heldOutFold.*NA, 1L or 2L")
  expect_error(initOnly(list(heldOutFold = "a")))
})

test_that("init schedules each event once", {
  for (p in list(list(), list(doFit = TRUE), list(mode = c("fit", "validate")), list(mode = "debug"))) {
    ev <- queued(initOnly(p))
    expect_false(anyDuplicated(ev) > 0L)
  }
})
