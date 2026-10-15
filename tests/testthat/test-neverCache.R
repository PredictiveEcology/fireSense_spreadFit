## `.neverCache` (read by SpaDES.core >= 3.2.1.9037) wins over every form of `.useCache`, including a
## `.globals$.useCache` that lists the event. `run` fits the model and writes the shared ledger through
## `CacheGeo()`; a cached `run` would replay an older fit without writing the ledger.
test_that(".neverCache is declared and SpaDES.core refuses to cache those events whatever .useCache says", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  declared <- md$parameters$default[[which(md$parameters$paramName == ".neverCache")]]
  expect_setequal(declared, c("run"))
  eventIsCached <- utils::getFromNamespace(".eventIsCached", "SpaDES.core")
  for (useCache in list(c("init", declared), TRUE)) {
    for (ev in declared)
      expect_false(eventIsCached(list(.useCache = useCache, .neverCache = declared), ev)$cache)
  }
  ## and the helper still caches what `.useCache` names that `.neverCache` does not list
  expect_true(eventIsCached(list(.useCache = "init", .neverCache = declared), "init")$cache)
})
