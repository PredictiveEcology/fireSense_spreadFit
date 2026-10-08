## 2026-09-29: fitSpread() called termsInDEoptim(), which named every non-formula parameter
## "logit1", "logit2", ... (maxAsymptote, inflectionPoint1, ...). fireSenseUtils::runDEoptim()
## (>= 0.2.3.9066) prints the real names from names(lower), so the module must not print its own.
test_that("fitSpread() does not call termsInDEoptim()", {
  expect_false("termsInDEoptim" %in% all.names(body(fitSpread)))
})
