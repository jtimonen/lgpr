library(lgpr)

# -------------------------------------------------------------------------

context("Methods for lgpmodel objects")

sim <- simulate_data(
  N = 4,
  t_data = seq(6, 36, by = 6),
  covariates = c(0, 2),
  t_observed = "after_1"
)
model <- create_model(y ~ gp(age) + zerosum(z) * gp(age), data = sim@data)

test_that("lgpmodel getters work", {
  n <- get_num_obs(model)
  df <- get_component_info(model)
  nams <- get_covariate_names(model)
  expect_equal(n, 24)
  expect_equal(dim(df), c(2,9))
  expect_equal(nams, "age, z")
})

test_that("prior summary works", {
  ps <- prior_summary(model, digits = 4)
  pars <- c("alpha[1]", "alpha[2]", "ell[1]", "ell[2]", "sigma[1]")
  expect_equal(ps$Parameter, pars)
})

test_that("model summary prints output", {
  expect_output(model_summary(model))
})

test_that("print_stan_input prints output", {
  expect_output(print_stan_input(model))
})
