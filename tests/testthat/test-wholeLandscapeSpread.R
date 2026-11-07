## wholeLandscapeSpread: one switch, passed to fireSenseUtils as `wholeLandscape` on every path that evaluates the
## objective; FALSE leaves the objective's arguments and the fit's cache key as they were.
## Helpers: fitArgs(), debugArgs() (helper-debugEvent.R), toySim() (helper-toyInputs.R).

wholeParams <- list(wholeLandscapeSpread = TRUE, spreadLink = "powerLaw", powerLawGamma = 1.9,
                    penaliseRunaways = FALSE)

test_that("FALSE passes the objective nothing new", {
  expect_false("wholeLandscape" %in% names(crpsObjectiveArgs(list(wholeLandscapeSpread = FALSE))))
  expect_false("wholeLandscape" %in% names(fitArgs(list(wholeLandscapeSpread = FALSE))))
})

test_that("TRUE reaches the fit (so the folds and the emulator, which take its arguments) and the debug event", {
  expect_true(fitArgs(wholeParams)$wholeLandscape)
  expect_true("wholeLandscape" %in% names(emulatorObjectiveArgs(c(fitArgs(wholeParams), list(lower = c(a = 1))))))
  expect_true(debugArgs(wholeParams)$wholeLandscape)
  expect_false(debugArgs()$wholeLandscape)
})

test_that("TRUE needs the powerLaw link and penaliseRunaways = FALSE", {
  expect_error(checkWholeLandscapeSpread(list(wholeLandscapeSpread = TRUE, spreadLink = "ceiling")), "powerLaw")
  expect_error(checkWholeLandscapeSpread(list(wholeLandscapeSpread = TRUE, spreadLink = "powerLaw",
                                              penaliseRunaways = TRUE)), "penaliseRunaways")
  expect_null(checkWholeLandscapeSpread(list(wholeLandscapeSpread = FALSE, spreadLink = "ceiling")))
  expect_error(fitArgs(modifyList(wholeParams, list(penaliseRunaways = TRUE))), "penaliseRunaways")
})
