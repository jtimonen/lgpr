library(lgpr)

# -------------------------------------------------------------------------

context("Prior sampling with different modeling options")

test_that("a model with uncertain disease age needs prior specified", {
  formula <- y ~ gp(age) + unc(id) * gp_vm(dis_age)
  data <- testdata_001
  reason <- "you must specify 'effect_time_info' in"
  expect_error(create_model(formula = formula, data = data), reason)
})

test_that("model with uncertain disease age can be created and fit", {
  formula <- y ~ gp(age) + unc(id) * gp_vm(dis_age)
  data <- testdata_001
  et <- list(backwards = FALSE, lower = 15, upper = 30, zero = 0)
  prior <- list(effect_time_info = et)
  model <- create_model(
    formula = formula,
    data = data,
    prior = prior,
    prior_only = TRUE
  )
  si <- get_stan_input(model)
  expect_equal(si$num_bt, 2)
  expect_equal(si$num_vm, 1)
  expect_equal(si$num_ns, 1)
  expect_equal(si$num_heter, 0)
  fit <- sample_model(
    model,
    iter = 2000,
    chains = 1,
    refresh = 0,
    seed = 123
  )
  expect_s4_class(fit, "lgpfit")
  p1 <- plot_warp(fit)
  p2 <- plot_effect_times(fit)
  p3 <- plot_pred(fit, fit_alpha = 0.02)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("model with heterogeneous disease effect can be created and fit", {
  formula <- y ~ gp(age) + het(id) * gp_ns(dis_age)
  data <- testdata_001
  model <- create_model(
    formula = formula,
    data = data,
    prior_only = TRUE
  )
  si <- get_stan_input(model)
  expect_equal(si$num_bt, 2)
  expect_equal(si$num_vm, 0)
  expect_equal(si$num_ns, 1)
  expect_equal(si$num_heter, 1)
  fit <- sample_model(
    model,
    iter = 2000,
    chains = 1,
    refresh = 0,
    seed = 123
  )
  expect_s4_class(fit, "lgpfit")
  p1 <- plot_beta(fit)
  expect_s3_class(p1, "ggplot")
})
