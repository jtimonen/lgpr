library(lgpr)

test_that("repvec works correctly", {
  expect_equal(
    repvec(3, 1),
    matrix(3)
  )
  expect_equal(
    repvec(3, 3),
    matrix(c(3, 3, 3), ncol = 1, nrow = 3)
  )
  expect_equal(
    repvec(c(2, 3, 4), 2),
    matrix(c(2, 3, 4, 2, 3, 4), ncol = 3, nrow = 2, byrow = TRUE)
  )
})
