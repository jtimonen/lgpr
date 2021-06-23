set.seed(3215)
library(lgpr)

# -------------------------------------------------------------------------

context("Base kernel functions")

test_that("base kernels work", {
  expect_equal(
    kernel_zerosum(c(1, 2), c(3, 2, 1), M = 3),
    matrix(c(-0.5, -0.5, 1.0, -0.5, 1.0, -0.5),
      nrow = 2, ncol = 3, byrow = TRUE
    )
  )
  expect_equal(
    kernel_bin(c(1, 2), c(3, 2, 1), pos_class = 2),
    matrix(c(0, 0, 0, 1, 0, 0), nrow = 2, ncol = 3, byrow = FALSE)
  )
  expect_equal(
    kernel_eq(-2, -2, ell = 20),
    matrix(1)
  )
  expect_equal(
    dim(kernel_ns(c(1, 1, 2), c(0, 1), ell = 1, a = 1)),
    c(3, 2)
  )
})

test_that("kernel_beta works correctly", {
  K <- kernel_beta(c(0.1, 0.5, 1.0), c(1, 2, 4), c(1, 1, 2, 2, 3, 3))
  expect_equal(dim(K), c(3, 6))
  expect_equal(K[1, 1], 0)
  expect_equal(K[2, 4], sqrt(0.1 * 0.1))
  expect_equal(K[3, 6], sqrt(0.5 * 1.0))
})

test_that("base kernels give errors when supposed to", {
  expect_error(
    kernel_eq(0, c(-1, 0, 1), ell = 0),
    "<ell> must be positive"
  )
  expect_error(
    kernel_zerosum(0, c(-1, 0, 1)),
    "is missing, with no default"
  )
  expect_error(
    kernel_ns(0, c(-1, 0, 1), ell = 1, alpha = -1, a = 1),
    "<alpha> must be non-negative"
  )
  expect_error(
    kernel_ns(0, c(-1, 0, 1), a = 1), # ell missing
    "is missing, with no default"
  )
})


# -------------------------------------------------------------------------

context("Compare R and Stan versions of kernel functions")
n1 <- 25
n2 <- 14
alpha <- exp(rnorm(1))
ell <- exp(rnorm(1))
num_categ <- 3

test_that("kernel_eq works similarly in R and Stan code", {
  x1 <- sort(rnorm(n = n1))
  x2 <- sort(rnorm(n = n2))
  K1 <- kernel_eq(x1, x2, alpha, ell)
  K2 <- STAN_kernel_eq(x1, x2, alpha, ell, get_stream())
  expect_equal(K1, K2)
  expect_equal(dim(K1), c(n1, n2))
})

test_that("kernel_cat works similarly in R and Stan code", {
  x1 <- sample.int(n = num_categ, size = n1, replace = TRUE)
  x2 <- sample.int(n = num_categ, size = n2, replace = TRUE)
  K1 <- kernel_cat(x1, x2)
  K2 <- STAN_kernel_cat(x1, x2, get_stream())
  expect_equal(K1, K2)
  expect_equal(dim(K1), c(n1, n2))
})

test_that("kernel_bin works similarly in R and Stan code", {
  x1 <- sample.int(n = num_categ, size = n1, replace = TRUE) - 1
  x2 <- sample.int(n = num_categ, size = n2, replace = TRUE) - 1
  K1 <- kernel_bin(x1, x2)
  K2 <- STAN_kernel_bin(x1, x2, get_stream())
  expect_equal(K1, K2)
  expect_equal(dim(K1), c(n1, n2))
})

test_that("kernel_zerosum works similarly in R and Stan code", {
  x1 <- sample.int(n = num_categ, size = n1, replace = TRUE)
  x2 <- sample.int(n = num_categ, size = n2, replace = TRUE)
  K1 <- kernel_zerosum(x1, x2, num_categ)
  K2 <- STAN_kernel_zerosum(x1, x2, num_categ, get_stream())
  expect_equal(K1, K2)
  expect_equal(dim(K1), c(n1, n2))
})

test_that("kernel_beta works similarly in R and Stan code", {
  beta <- c(0.3, 0.1, 0.5, 1.0, 0.3)
  n_cases <- length(beta)
  idx1_expand <- sample.int(n = n_cases + 1, size = n1, replace = TRUE)
  idx2_expand <- sample.int(n = n_cases + 1, size = n2, replace = TRUE)
  # idx = 1 means control, 2-6 are cases with one beta param each
  K1 <- kernel_beta(beta, idx1_expand, idx2_expand)
  K2 <- STAN_kernel_beta(beta, idx1_expand, idx2_expand, get_stream())
  expect_equal(K1, K2)
  expect_equal(dim(K1), c(n1, n2))
})

test_that("kernel_varmask works similarly in R and Stan code", {
  x1 <- sort(rnorm(n = n1))
  x2 <- sort(rnorm(n = n2))
  a <- 0.3 + runif(1)
  vm_params <- runif(n = 2)
  K1 <- kernel_varmask(x1, x2, a, vm_params)
  K2 <- STAN_kernel_varmask(x1, x2, a, vm_params, get_stream())
  expect_equal(K1, K2)
  expect_equal(dim(K1), c(n1, n2))
})

# -------------------------------------------------------------------------

context("Compare diagonal and full versions of Stan kernel functions")
n1 <- 27
alpha <- exp(rnorm(1))
ell <- exp(rnorm(1))
num_categ <- 3

test_that("STAN_kernel_eq_diag works correctly", {
  x <- sort(rnorm(n = n1))
  K1 <- STAN_kernel_eq(x, x, alpha, ell, get_stream())
  D1 <- diag(K1)
  D2 <- STAN_kernel_eq_diag(length(x), alpha, get_stream())
  expect_equal(D1, D2)
  expect_equal(length(D1), n1)
})

test_that("STAN_kernel_const_diag works correctly (zerosum kernel)", {
  ktype <- 0 # 0 = zs, 1 = cat, 2 = bin
  x <- sample.int(n = num_categ, size = n1, replace = TRUE)
  K1 <- STAN_kernel_const(x, x, ktype, 0, get_stream())
  D1 <- diag(K1)
  D2 <- STAN_kernel_const_diag(x, ktype, get_stream())
  expect_equal(D1, D2)
  expect_equal(length(D1), n1)
  expect_equal(prod(D1), 1.0) # should be all ones
})

test_that("STAN_kernel_const_diag works correctly (categorical kernel)", {
  ktype <- 1 # 0 = zs, 1 = cat, 2 = bin
  x <- sample.int(n = num_categ, size = n1, replace = TRUE)
  K1 <- STAN_kernel_const(x, x, ktype, num_categ, get_stream())
  D1 <- diag(K1)
  D2 <- STAN_kernel_const_diag(x, ktype, get_stream())
  expect_equal(D1, D2)
  expect_equal(length(D1), n1)
  expect_equal(prod(D1), 1.0) # should be all ones
})

test_that("STAN_kernel_const_diag works correctly (binary mask kernel)", {
  ktype <- 2 # 0 = zs, 1 = cat, 2 = bin
  x <- sample.int(n = num_categ, size = n1, replace = TRUE) - 1
  K1 <- STAN_kernel_const(x, x, ktype, 0, get_stream())
  D1 <- diag(K1)
  D2 <- STAN_kernel_const_diag(x, ktype, get_stream())
  expect_equal(D1, D2) # should have ones and zeros
  expect_equal(length(D1), n1)
})

test_that("STAN_kernel_beta_diag works correctly", {
  beta <- c(0.29, 0.1, 0.59, 1.0, 0.43)
  n_cases <- length(beta)
  idx_expand <- sample.int(n = n_cases + 1, size = n1, replace = TRUE)
  # idx = 1 means control, 2-6 are cases with one beta param each
  K1 <- kernel_beta(beta, idx_expand, idx_expand)
  D1 <- diag(K1)
  D2 <- STAN_kernel_beta_diag(beta, idx_expand, get_stream())
  expect_equal(D1, D2)
  expect_equal(length(D1), n1)
})

test_that("kernel_varmask works similarly in R and Stan code", {
  x <- sort(rnorm(n = n1))
  a <- 0.23 + runif(1)
  vm_params <- runif(n = 2)
  K1 <- kernel_varmask(x, x, a, vm_params)
  D1 <- diag(K1)
  D2 <- STAN_kernel_varmask_diag(x, a, vm_params, get_stream())
  expect_equal(D1, D2)
  expect_equal(length(D1), n1)
})
