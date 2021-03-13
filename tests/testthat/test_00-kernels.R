library(lgpr)

# -------------------------------------------------------------------------

context("Kernel functions")

test_that("base kernels work", {
  expect_equal(
    kernel_zerosum(c(1, 2), c(3, 2, 1), M = 3, alpha = 1),
    matrix(c(-0.5, -0.5, 1.0, -0.5, 1.0, -0.5),
           nrow = 2, ncol = 3, byrow = TRUE
    )
  )
  expect_equal(
    kernel_bin(c(1, 2), c(3, 2, 1), pos_class = 2),
    matrix(c(0, 0, 0, 1, 0, 0), nrow = 2, ncol = 3, byrow = FALSE)
  )
  expect_equal(
    kernel_se(-2, -2, ell = 20),
    matrix(1)
  )
  expect_equal(
    dim(kernel_ns(c(1, 1, 2), c(0, 1), ell = 1, a = 1)),
    c(3, 2)
  )
})

test_that("kernel_beta works correctly", {
  K <- kernel_beta(c(0.1, 0.5, 1.0), c(1, 2, 3), c(1, 1, 2, 2, 3, 3))
  expect_equal(dim(K), c(3, 6))
  expect_equal(K[1, 1], 0.1)
  expect_equal(K[2, 4], 0.5)
  expect_equal(K[3, 6], 1.0)
})

test_that("base kernels give errors when supposed to", {
  expect_error(
    kernel_se(0, c(-1, 0, 1), ell = 0), 
    "<ell> must be positive"
  )
  expect_error(
    kernel_zerosum(0, c(-1, 0, 1), alpha = 1),
    "is missing, with no default"
  )
  expect_error(
    kernel_ns(0, c(-1, 0, 1), ell = 1, alpha = -1, a = 1),
    "<alpha> must be non-negative"
  )
  expect_error(
    kernel_ns(0, c(-1, 0, 1), alpha = 1, a = 1),
    "is missing, with no default"
  )
})
