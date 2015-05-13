context("Utilities")

test_that("deduplicate removes duplicates", {
  a <- list(a = letters[1:2], b = letters[1:2], c = letters[2:3])
  expect_identical(deduplicate(a), a[c(1, 3)])
}) 

