library(lgpr)
set.seed(123)

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
  expect_equal(dim(df), c(2, 9))
  expect_equal(nams, "age, z")
})

test_that("prior summary works", {
  ps <- prior_summary(model, digits = 4)
  pars <- c("alpha[1]", "alpha[2]", "ell[1]", "ell[2]", "(sigma[1])^2")
  expect_equal(ps$Parameter, pars)
})

test_that("model summary prints output", {
  expect_output(model_summary(model))
  expect_output(show(model))
  expect_output(show(model@model_formula))
})

test_that("print_stan_input prints output", {
  expect_output(print_stan_input(model))
})


# -------------------------------------------------------------------------

context("Methods for lgpsim objects")

test_that("simulated data can be plotted", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(1, 2)
  )
  p <- sim_plot(dat, i_test = c(1, 2, 3), ncol = 4)
  expect_s3_class(p, "ggplot")
})

test_that("simulated data with disease effect can be plotted", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  p <- sim_plot(dat, i_test = c(1, 2, 3), ncol = 4)
  expect_s3_class(p, "ggplot")
})

test_that("show method for simualated data prints output", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  expect_output(show(dat))
})


# -------------------------------------------------------------------------

context("Methods for lgpfit objects")

sim <- simulate_data(
  N = 4,
  t_data = seq(6, 36, by = 6),
  covariates = c(0, 2),
  t_observed = "after_1"
)
data <- sim@data
data$y <- data$y + 5

test_that("model fit can be visualized with data on original scale", {
  suppressWarnings({
    fit <- lgp(y ~ gp(age) + zerosum(z) * gp(age) + gp_warp_vm(diseaseAge),
      data = data, chains = 1, iter = 300, refresh = 0
    )
    p1 <- plot_fit(fit, data)
    p2 <- plot_fit(fit, data, draws = c(2, 3))
    p3 <- plot_fit(fit, data, draws = 3)
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    expect_s3_class(p3, "ggplot")
  })
})

test_that("fit summary prints output", {
  suppressWarnings({
    fit <- lgp(y ~ gp(age) + zerosum(z) * gp(age),
      data = data, chains = 1, iter = 100, refresh = 0
    )
    expect_output(fit_summary(fit))
    expect_output(show(fit))
  })
})
