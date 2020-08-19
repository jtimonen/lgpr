library(lgpr)

context("Stan stream")
STREAM <- get_stream()

# -------------------------------------------------------------------------

context("Stan base kernels: zero-sum kernel")

test_that("zero-sum kernel works correctly", {
  M <- 2
  x <- c(1, 1, 2)
  a <- c(
    1, 1, -1,
    1, 1, -1,
    -1, -1, 1
  )
  K_expect <- matrix(a, 3, 3, byrow = TRUE)
  K <- STAN_kernel_base_zerosum(x, x, M, STREAM)
  expect_equal(K, K_expect)
})

test_that("zero-sum kernel works similarly as reference", {
  M <- 3
  x <- sample.int(M, size = 8, replace = TRUE)
  expect_equal(
    STAN_kernel_base_zerosum(x, x, M, STREAM),
    sim_kernel_zerosum(x, x, M)
  )
})

test_that("zero-sum kernel errors if number of categories is one", {
  M <- 1
  x <- sample.int(M, size = 8, replace = TRUE)
  expect_error(STAN_kernel_base_zerosum(x, x, M, STREAM))
})


context("Stan base kernels: ordinary categorical kernel")

test_that("categorical kernel works correctly", {
  x <- c(1, 1, 2)
  a <- c(
    1, 1, 0,
    1, 1, 0,
    0, 0, 1
  )
  K_expect <- matrix(a, 3, 3, byrow = TRUE)
  K <- STAN_kernel_base_cat(x, x, STREAM)
  expect_equal(K, K_expect)
})


context("Stan base kernels: binary mask kernel")

test_that("binary mask kernel works correctly", {
  x <- c(0, 0, 1)
  a <- c(
    1, 1, 0,
    1, 1, 0,
    0, 0, 0
  )
  K_expect <- matrix(a, 3, 3, byrow = TRUE)
  K <- STAN_kernel_base_bin_mask(x, x, STREAM)
  expect_equal(K, K_expect)
})


context("Stan base kernels: variance mask kernel")

test_that("variance mask kernel works correctly", {
  x <- c(12, 0, 12)
  stp <- 0.2
  vm_params <- c(0.05, 0.6)
  v <- c(
    0.9755191, 0.9382995, 0.9755191,
    0.9382995, 0.9025000, 0.9382995,
    0.9755191, 0.9382995, 0.9755191
  )
  K_expect <- matrix(v, 3, 3, byrow = TRUE)
  K <- STAN_kernel_base_var_mask(x, x, stp, vm_params, STREAM)
  expect_equal(K, K_expect)
})


context("Stan base kernels: variance mask kernel")

test_that("variance mask kernel works similarly as reference", {
  x <- c(-24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12)
  stp <- 1.0
  vm_params <- c(0.05, 0.6)
  K_stan <- STAN_kernel_base_var_mask(x, x, stp, vm_params, STREAM)
  K_r <- sim_kernel_var_mask(x, x, vm_params, stp)
  expect_equal(K_stan, K_r)
})

test_that("variance mask kernel errors if steepness is not valid", {
  x <- c(-24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12)
  vm_params <- c(0.05, 0.6)
  expect_error(STAN_kernel_base_var_mask(x, x, 0, vm_params, STREAM))
  expect_error(STAN_kernel_base_var_mask(x, x, -1.0, vm_params, STREAM))
  expect_error(STAN_kernel_base_var_mask(x, x, NaN, vm_params, STREAM))
  expect_error(STAN_kernel_base_var_mask(x, x, Inf, vm_params, STREAM))
})

test_that("variance mask kernel errors if <vm_params> is not valid", {
  x <- c(-24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12)
  vm_params <- c(0.05, 0.6)
  expect_error(STAN_kernel_base_var_mask(x, x, 1, c(-0.05, 0.6), STREAM))
  expect_error(STAN_kernel_base_var_mask(x, x, 1, c(1.6, 0.6), STREAM))
  expect_error(STAN_kernel_base_var_mask(x, x, 1, c(0.05, NaN), STREAM))
  expect_error(STAN_kernel_base_var_mask(x, x, 1, c(NaN, 0.6), STREAM))
})
