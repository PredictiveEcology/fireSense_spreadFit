## ELF 4.2.1 on the GPSC: the ledger read before the fit had its neighbours' rows but not its own.
## fireSense_dataPrepFit built the per-ELF species and land-cover lists from those rows; the fit then wrote
## 4.2.1's row, and prediction in the same run (fireSense_dataPrepPredict's ELFfuelSets()) stopped with
## "sppEquivs has no entry for ELF(s) 4.2.1 of sim$studyAreaWithSpreadParams". The `run` event now
## rebuilds the lists from the ledger it returns.

neighbourSppEquiv <- function() data.table::data.table(LandR = c("Abie_bal", "Pice_gla"), FuelClass = c("class1", "class2"))

## CacheGeo(action = "update") returns the ledger merged with the new row: here, a neighbour and the new row.
mockLedgerWithNeighbour <- function(sim, rec = new.env()) {
  mockFitAndLedger(sim, rec)
  mockInModule(sim, CacheGeo = function(cloudFolderID, targetFile, domain, destinationPath, FUN, ..., action) {
    rec$geoArgs <- list(action = action)
    own <- list(...)$studyAreaFireSense
    nb <- own
    nb$polygonID <- "9.8"
    nb$sppEquiv <- list(neighbourSppEquiv())
    nb$nonForestedLCCGroups <- list(list(nf_nb = 50))
    nb$missingLCCgroup <- "nf_nb"
    rbind(nb, own)
  })
  rec
}

## What fireSense_dataPrepFit hands over when the ledger has only the neighbour: its row and its lists.
neighbourOnlyObjects <- function() {
  list(studyAreaWithSpreadParams = toyLedger("9.8"),
       sppEquivs = list(`9.8` = neighbourSppEquiv()), sppColorVects = list(`9.8` = c(Abie_bal = "#000001")),
       sppNameVectors = list(`9.8` = c("Abie_bal", "Pice_gla")),
       nonForestedLCCGroupsList = list(`9.8` = list(nf_nb = 50)), missingLCCgroupList = list(`9.8` = "nf_nb"))
}

ledgerLists <- c("sppEquivs", "sppColorVects", "sppNameVectors", "nonForestedLCCGroupsList", "missingLCCgroupList")

test_that("after a fit whose ELF had only neighbour rows, the per-ELF lists cover every ledger row", {
  sim <- toySim(objects = neighbourOnlyObjects())
  mockLedgerWithNeighbour(sim)
  sim <- suppressMessages(SpaDES.core::spades(sim))
  ids <- as.character(sim$studyAreaWithSpreadParams$polygonID)
  expect_identical(ids, c("9.8", "9.9"))
  for (l in ledgerLists) expect_identical(names(sim[[l]]), ids, label = l)
  ## as fireSense_dataPrepPredict's ELFfuelSets() finds each row's entry
  expect_false(anyNA(match(ids, names(sim$sppEquivs))))
  expect_identical(sim$sppEquivs[["9.9"]]$LandR, sim$sppEquiv$LandR)
  expect_identical(sim$nonForestedLCCGroupsList[["9.9"]], sim$nonForestedLCCGroups)
  expect_identical(sim$missingLCCgroupList[["9.9"]], sim$missingLCCgroup)
  expect_identical(sim$sppEquivs[["9.8"]]$LandR, c("Abie_bal", "Pice_gla"))
})

test_that("the rebuilt lists are outputs, so a cached `run` event restores them", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_true(all(ledgerLists %in% md$outputObjects$objectName))
  withr::local_options(reproducible.cachePath = toyPaths()$cachePath)
  first <- toySim(params = list(.useCache = "run"), objects = neighbourOnlyObjects())
  mockLedgerWithNeighbour(first)
  first <- suppressMessages(SpaDES.core::spades(first))
  second <- toySim(params = list(.useCache = "run"), objects = neighbourOnlyObjects())
  rec <- mockLedgerWithNeighbour(second)
  second <- suppressMessages(SpaDES.core::spades(second))
  expect_null(rec$geoArgs)                                   # the event came from the cache
  for (l in ledgerLists) expect_identical(names(second[[l]]), c("9.8", "9.9"), label = l)
})
