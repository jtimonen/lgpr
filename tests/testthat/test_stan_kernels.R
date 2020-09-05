library(lgpr)
STREAM <- get_stream()

# -------------------------------------------------------------------------

context("Stan base kernels")

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



# -------------------------------------------------------------------------

context("Stan high-level functions")

# Create test input
sim <- simulate_data(
  N = 4,
  t_data = seq(6, 36, by = 6),
  covariates = c(0, 1, 2, 3),
  lengthscales = rep(12, 5),
  relevances = rep(1, 6),
  t_jitter = 0.5
)

# Model
m <- create_model(y ~ zerosum(id) * gp(age) + gp_warp_vm(diseaseAge) +
  categ(z) + gp(age) + gp(x),
data = sim@data
)

# Input
input <- m@stan_input
x_cat <- input$x_cat
x_mask <- input$x_cont_mask
num_levels <- input$x_cat_num_levels
comp <- input$components
K_const <- kernel_const_all(x_cat, x_cat, x_mask, x_mask, num_levels, comp)

# Params
alpha <- c(1, 1, 1, 1, 1)
ell <- c(1, 1, 1, 1)
x <- input$x_cont
x_unnorm <- input$x_cont_unnorm
vm_params <- input$vm_params
ix <- input$idx_expand

# All kernels
teff_zero <- dollar(input, "teff_zero")
K <- kernel_all(
  K_const, comp, x, x, x_unnorm, x_unnorm,
  alpha, ell, 0.5, list(), list(),
  vm_params, ix, ix, teff_zero
)

test_that("kernel_const_all works correctly", {
  expect_equal(length(K_const), 5)
  expect_equal(dim(K_const[[1]]), c(24, 24))
  S1 <- sum(K_const[[1]])
  expect_lt(abs(S1), 1e-6)
  expect_equal(sum(K_const[[4]]), 24 * 24)
  expect_equal(sum(K_const[[5]]), 24 * 24)
})

test_that("kernel_all works correctly", {
  expect_equal(length(K), 5)
  expect_equal(dim(K[[1]]), c(24, 24))
})

test_that("gp_posterior works correctly", {
  y <- sim@data$y
  fp <- gp_posterior(K, y, 1e-6, 1.0)
  expect_equal(length(fp), 12)

  # test that componentwise means sum to total mean
  f_sum <- fp[[1]]
  for (j in 2:5) {
    f_sum <- f_sum + fp[[j]]
  }
  diff <- f_sum - fp[[6]]
  expect_lt(max(abs(diff)), 1e-6)
})
