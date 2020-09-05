library(lgpr)

# -------------------------------------------------------------------------

context("Methods for lgpmodel objects")

et <- list(lower = 10, zero = 0, backwards = FALSE, upper = 20)
model <- create_model(y ~ uncrt(id) * heter(id) * gp_warp(dis_age) +
  zerosum(sex) * gp(age),
prior = list(effect_time_info = et),
data = testdata_001
)

test_that("lgpmodel getters work", {
  n <- get_num_obs(model)
  df <- get_component_info(model)
  nams <- get_covariate_names(model)
  expect_equal(n, 24)
  expect_equal(dim(df), c(2, 9))
  expect_equal(nams, "dis_age, age, id, sex")
})

test_that("prior summary works", {
  ps <- prior_summary(model, digits = 4)
  expect_equal(length(ps$Parameter), 9)
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
set.seed(123)

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
