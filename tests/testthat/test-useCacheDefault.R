## A project that sets no `.useCache` gets the module's default. It must cache what the
## project script used to switch on for every module: `.inputObjects`, `init`, `estimateThreshold` and
## `spreadFitPrepare`. (`initPlot` and `checkData` of that list are not events or cache switches here.)

md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
useCacheDefault <- md$parameters$default[[which(md$parameters$paramName == ".useCache")]]
code <- readLines(file.path(moduleRoot, paste0(moduleName, ".R")))

test_that("the default .useCache is every cached event of this module", {
  expect_setequal(useCacheDefault, c(".inputObjects", "init", "estimateThreshold", "spreadFitPrepare"))
})

test_that("each name in the default .useCache is an event or .inputObjects of this module", {
  eventsScheduled <- unique(sub('.*moduleName, "([[:alnum:]]+)".*', "\\1", grep('scheduleEvent\\(.*moduleName, "', code, value = TRUE)))
  expect_true(any(grepl("^\\.inputObjects <- function", code)))
  expect_true(all(setdiff(useCacheDefault, c("init", ".inputObjects")) %in% eventsScheduled))
})

test_that("the module does not read the phase parameter (fireSenseUtils::runJob() does)", {
  allCode <- c(code, unlist(lapply(list.files(file.path(moduleRoot, "R"), full.names = TRUE), readLines)))
  expect_false(any(grepl("(P\\(sim\\)|Par|params\\(sim\\)\\[\\[[^]]+\\]\\])\\$phase\\b|\\[\\[\"phase\"\\]\\]", allCode)))
})
