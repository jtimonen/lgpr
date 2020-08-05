library(lgpr)
library(rstan)
stanmodel <- lgpr::get_stan_model()
rstan::expose_stan_functions(stanmodel)

# Create data for testing
N <- 4
set.seed(123)
data <- lgpr::simulate_data(
  N = N,
  t_data = c(6, 12, 18, 24, 30),
  covariates = c(0, 1, 2, 2)
)$data
num_obs <- dim(data)[1]
x <- list(
  data[, 1], data[, 2], data[, 3],
  data[, 4], data[, 5], data[, 6], rep(0, num_obs)
)

x_num_cat <- c(4, 0, 0, 0, 2, 2)
x_caseid <- x[[1]] - 2
x_caseid[which(x_caseid < 0)] <- 0
ctypes <- c(0, 1, 2, 2, 3)
ktypes <- c(0, 0, 0, 0, 1)
num_comps <- length(ctypes)

cov1 <- c(1, num_comps, 5, 6, num_comps)
cov2 <- c(0, 2, 2, 2, 3)
components <- list(ctypes, ktypes, cov1, cov2)
num_cases <- sum(x_caseid > 0)

# Run tests
context("STAN_kernel_fixed_all")

test_that("STAN_kernel_fixed_all returns the correct number of components", {
  KF <- STAN_kernel_fixed_all(x, x_num_cat, x_caseid, components)
  expect_equal(length(KF), num_comps)
})

test_that("STAN_kernel_fixed_all returns matrices of correct size", {
  KF <- STAN_kernel_fixed_all(x, x_num_cat, x_caseid, components)
  K1 <- KF[[1]]
  expect_equal(dim(K1)[1], num_obs)
  expect_equal(dim(K1)[2], num_obs)
})

test_that("STAN_kernel_fixed_all returns matrices with correct values", {
  KF <- STAN_kernel_fixed_all(x, x_num_cat, x_caseid, components)
  expect_equal(KF[[1]][1, 1], 1)
  expect_equal(KF[[2]][1, 1], 0)
})

context("STAN_count_num_params")

test_that("STAN_count_num_params works", {
  num_params <- STAN_count_num_params(components, 1, num_cases)
  expect_equal(length(num_params), 6)
})
