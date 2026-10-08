## hasPreRunFitForThisPolygon(): "is there a ledger row for THIS polygon", not "is there a ledger".
## It reads three elements with `$`, so a plain list stands in for the simList.

led <- function(ids) toyLedger(ids)

test_that("no ledger, or something that is not a table, is not a fit", {
  expect_false(hasPreRunFitForThisPolygon(list(.ELFind = "4.3")))
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = "4.3", .ELFind = "4.3")))
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = list(polygonID = "4.3"),
                                               .ELFind = "4.3")))
})

test_that("an empty ledger, or one without polygonID, is not a fit", {
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = led("4.3")[0, ],
                                               .ELFind = "4.3")))
  noKey <- led("4.3"); names(noKey)[names(noKey) == "polygonID"] <- "id"
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = noKey, .ELFind = "4.3")))
})

test_that("only a row for this polygon counts: neighbours' rows do not", {
  neighbours <- led(c("4.1", "4.2", "14.3"))
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = neighbours, .ELFind = "4.3")))
  expect_true(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = neighbours, .ELFind = "14.3")))
  ## whole-string match: "4.3" is not found in "14.3" and "4" is not found in "4.1"
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = neighbours, .ELFind = "4")))
})

test_that("a plain data.frame ledger and a numeric polygonID both work", {
  df <- data.frame(polygonID = c(4.1, 4.3))
  expect_true(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = df, .ELFind = "4.3")))
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = df, .ELFind = "4.2")))
})

test_that(".ELFind is the key; .runName is used only when .ELFind is absent", {
  l <- led("4.3")
  expect_true(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = l, .runName = "4.3")))
  expect_true(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = l, .ELFind = character(0),
                                              .runName = "4.3")))
  ## .ELFind wins over a matching .runName
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = l, .ELFind = "9.9",
                                               .runName = "4.3")))
  ## neither supplied
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = l)))
  ## only the first element of a vector id is used
  expect_true(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = l, .ELFind = c("4.3", "9.9"))))
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = l, .ELFind = c("9.9", "4.3"))))
})

test_that("an NA polygonID or NA id is never a match", {
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = data.frame(polygonID = NA),
                                               .ELFind = "4.3")))
  expect_false(hasPreRunFitForThisPolygon(list(studyAreaWithSpreadParams = led("4.3"),
                                               .ELFind = NA_character_)))
})
