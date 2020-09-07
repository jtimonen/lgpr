library(lgpr)

# -------------------------------------------------------------------------

context("Prior sampling")

test_that("a model with nb likelihood can be sampled", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "nb",
      data = dat,
      iter = 600,
      chains = 1,
      refresh = 0,
    )
    expect_s4_class(fit, "lgpfit")
    p1 <- plot_draws(fit)
    p2 <- plot_draws(fit, regex_pars = "f_latent")
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")

    a <- get_f(fit)
    expect_equal(dim(a$f$`gp(age)`), c(300, 16))
  })
})


test_that("a model with uncertain disease age needs prior specified", {
  formula <- y ~ gp(age) + uncrt(id) * gp_warp_vm(dis_age)
  data <- testdata_001
  reason <- "you must specify 'effect_time_info' in the prior list"
  expect_error(create_model(formula = formula, data = data), reason)
})

test_that("model with uncertain disease age can be created and fit", {
  formula <- y ~ gp(age) + uncrt(id) * gp_warp_vm(dis_age)
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
  p3 <- plot_fit(fit, data = data)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("model with heterogeneous disease effect can be created and fit", {
  formula <- y ~ gp(age) + heter(id) * gp_warp(dis_age)
  data <- testdata_001
  model <- create_model(
    formula = formula,
    data = data,
    prior = prior,
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
