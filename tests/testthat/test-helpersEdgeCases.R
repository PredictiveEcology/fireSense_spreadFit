## Edge cases of bestParamSets() not covered by the files that introduced them.

test_that("bestParamSets(): n sets how many members are recorded, best first", {
  skip_if_not("bestByReplicatedMean" %in% getNamespaceExports("fireSenseUtils"))
  DE <- toyDE(nPar = 2L)                       # member k has p1 = k and value 8 - k
  out <- bestParamSets(DE, c("p1", "p2"), n = 2L)
  expect_equal(out$params$p1, c(7, 6))
  expect_equal(out$objFunVal, c(1, 2))
  ## only the LAST block's population is read
  DE[[1]]$member$pop[, 1] <- 99
  expect_equal(bestParamSets(DE, c("p1", "p2"), n = 1L)$params$p1, 7)
})
