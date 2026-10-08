## .elfSeed(): the threshold seed is a pure function of the polygon id.

test_that(".elfSeed() is 1 + (sum of byte * position * 7919) mod 1e6", {
  ## "A" is byte 65: 65 * 1 * 7919 = 514735
  expect_identical(.elfSeed("A"), 514736L)
  ## "AB": 7919 * (65 * 1 + 66 * 2) = 7919 * 197 = 1560043; mod 1e6 = 560043
  expect_identical(.elfSeed("AB"), 560044L)
  ## position matters: "BA" = 7919 * (66 + 130) = 1552124
  expect_identical(.elfSeed("BA"), 552125L)
  ## "4.3": bytes 52, 46, 51 -> 7919 * (52 + 92 + 153) = 7919 * 297 = 2351943
  expect_identical(.elfSeed("4.3"), 351944L)
})

test_that(".elfSeed() uses the first 64 bytes only", {
  base64 <- strrep("a", 64)
  expect_identical(.elfSeed(paste0(base64, "x")), .elfSeed(paste0(base64, "y")))
  expect_false(identical(.elfSeed(strrep("a", 63)), .elfSeed(base64)))
})

test_that(".elfSeed() refuses anything but one non-empty string", {
  expect_error(.elfSeed(4.3), "is.character(elf) is not TRUE", fixed = TRUE)
  expect_error(.elfSeed(c("4.3", "4.4")), "length(elf) == 1L is not TRUE", fixed = TRUE)
  expect_error(.elfSeed(""), "nzchar(elf) is not TRUE", fixed = TRUE)
  expect_error(.elfSeed(NULL), "is.character(elf) is not TRUE", fixed = TRUE)
})
