
context("stan helper functions")

library(lgpr)
library(rstan)
stanmodel <- lgpr::get_stan_model()
rstan::expose_stan_functions(stanmodel)

test_that("input warping function works similarly in Stan and R", {
  a <- exp(stats::rnorm(1))
  x <- c(-2, 1, 0, 1, 2)
  expect_equal(
    STAN_warp_input(x, a),
    warp_input(x, a, 0, 1)
  )
})

test_that("variance masking function works similarly in Stan and R", {
  a <- 0.6
  x <- c(-5, 0, 5)
  expect_equal(
    STAN_var_mask(x = x, a = a),
    var_mask(x = x, a = a)
  )
})

test_that("STAN_get_x_tilde works properly", {
  x_disAge <- c(-24, -12, 0, 12, -24, -12, 0, 12, 0, 0, 0, 0, 0, -12, 0, 12, 16)
  T_effect <- c(-1, 2, 10)
  T_observed <- c(0, 6, 12)
  mapping <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 14, 15, 16, 17), ncol = 4, nrow = 3, byrow = TRUE)
  map <- list(mapping[1, ], mapping[2, ], mapping[3, ])
  lengths <- c(4, 4, 4)
  expected <- c(-23, -11, 1, 13, -20, -8, 4, 16, 0, 0, 0, 0, 0, -10, 2, 14, 18)
  expect_equal(
    STAN_get_x_tilde(x_disAge, T_effect, T_observed, map, lengths),
    expected
  )
})

context("stan kernel functions")

test_that("zero-sum kernel works similarly in R and Stan", {
  M <- 3
  x <- sample.int(M, size = 8, replace = TRUE)
  expect_equal(
    STAN_K_zerosum(x, x, M),
    kernel_zerosum(x, x, M)
  )
})

test_that("binary mask kernel works similarly in R and Stan", {
  x <- sample.int(2, size = 8, replace = TRUE) - 1
  expect_equal(
    STAN_K_bin(x, x, 1),
    kernel_bin(x, x)
  )
})

test_that("variance mask kernel works similarly in R and Stan", {
  x <- c(-24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12)
  stp <- 1.0
  vm_params <- c(0.05, 0.6)
  expect_equal(
    STAN_K_var_mask(x, stp, vm_params),
    kernel_var_mask(x, x, vm_params, stp)
  )
})
