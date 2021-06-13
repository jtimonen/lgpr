library(lgpr)
source("helpers/SW.R")

# -------------------------------------------------------------------------

context("Optimizing MAP parameters")

test_that("optimize_model can optimize MAP parameters", {
  dat <- testdata_001
  model <- create_model(y ~ gp(age), dat)
  fit <- optimize_model(model, iter = 10, seed = 123)
  expect_equal(class(fit), "list")
})


# -------------------------------------------------------------------------

N_ITER <- 23
N_CHAINS <- 1
SEED <- 992

context("Posterior sampling (f marginalized)")

DAT <- testdata_001

my_prior <- list(
  alpha = normal(1, 0.1),
  ell = normal(1, 0.1),
  sigma = normal(1, 0.1)
)


test_that("lgp() can do posterior sampling (f marginalized)", {
  SW({
    fit <- lgp(
      y ~ id + age,
      data = DAT,
      iter = N_ITER,
      prior = my_prior,
      chains = N_CHAINS,
      refresh = 0,
      seed = 123
    )
  })

  expect_s4_class(fit, "lgpfit")
  p1 <- plot(fit)
  p2 <- plot_draws(fit, type = "trace")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_error(plot_warp(fit), "model does not have input warping parameters")
  expect_error(plot_beta(fit), "there are no heterogeneous effects")
  expect_error(plot_effect_times(fit), "there are no uncertain effect times")
  expect_output(show(fit))
  expect_true(!is_f_sampled(fit))
})


test_that("verbose mode can be used in lgp()", {
  dat <- testdata_001
  expect_output(
    suppressWarnings({
      lgp(
        formula = y ~ gp(age),
        data = DAT,
        iter = N_ITER,
        chains = N_CHAINS,
      )
    })
  )
})

# -------------------------------------------------------------------------

context("Posterior sampling with different obs models (f latent)")

test_that("f can be sampled with gaussian likelihood", {
  SW({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex),
      sample_f = TRUE,
      data = DAT,
      iter = N_ITER,
      chains = N_CHAINS,
      refresh = 0,
    )
  })
  expect_s4_class(fit, "lgpfit")
  expect_true(is_f_sampled(fit))
})

# Positive integer data
NEWDAT <- DAT
NEWDAT$y <- round(exp(NEWDAT$y))

test_that("f can be sampled with poisson likelihood", {
  SW({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "poisson",
      data = NEWDAT,
      iter = N_ITER,
      chains = N_CHAINS,
      refresh = 0,
    )
  })
  expect_s4_class(fit, "lgpfit")
  expect_equal(get_obs_model(fit), "poisson")
  expect_true(is_f_sampled(fit))
})

test_that("f can be sampled with nb likelihood", {
  SW({
    fit <- lgp(
      formula = y ~ gp_ns(age) + categ(sex) * gp(age),
      likelihood = "nb",
      data = NEWDAT,
      iter = N_ITER,
      chains = N_CHAINS,
      refresh = 0,
    )
  })
  expect_s4_class(fit, "lgpfit")
  p1 <- plot(fit)
  p2 <- plot_draws(fit, regex_pars = "f_latent")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_equal(get_obs_model(fit), "nb")
  expect_true(is_f_sampled(fit))
})

test_that("f can be sampled with binomial likelihood", {
  SW({
    fit <- lgp(
      formula = y ~ gp(age) + zs(sex) * gp_vm(age),
      likelihood = "binomial",
      data = NEWDAT,
      iter = N_ITER,
      chains = N_CHAINS,
      refresh = 0,
      num_trials = 10
    )
  })
  expect_s4_class(fit, "lgpfit")
  expect_equal(get_obs_model(fit), "binomial")
  expect_true(is_f_sampled(fit))
})

test_that("f can be sampled with beta-binomial likelihood", {
  SW({
    fit <- lgp(
      formula = y ~ gp(age) + categ(id) * gp(age) + categ(id),
      likelihood = "bb",
      data = NEWDAT,
      iter = N_ITER,
      chains = N_CHAINS,
      refresh = 0,
      num_trials = 10
    )
  })
  ci <- get_component_encoding(fit)
  expect_equal(as.numeric(ci[, 1]), c(1, 2, 0)) # types
  expect_equal(as.numeric(ci[, 2]), c(0, 1, 1)) # kernels
  expect_s4_class(fit, "lgpfit")
  expect_output(show(fit@model))
  expect_equal(get_obs_model(fit), "bb")
  expect_true(is_f_sampled(fit))
})

# -------------------------------------------------------------------------

context("Prior sampling with prior_only = TRUE")

test_that("Model with uncertain effect time (f marginalized)", {
  formula <- y ~ gp(age) + unc(id) * gp_vm(dis_age)
  data <- testdata_001
  et <- list(backwards = FALSE, lower = 15, upper = 30, zero = 0)
  prior <- list(effect_time_info = et, wrp = igam(14, 5))

  # NOTE: Not suppressing warnings here!
  fit <- lgp(
    formula = formula,
    data = data,
    prior = prior,
    prior_only = TRUE,
    iter = 2000,
    chains = 1,
    refresh = 0,
    seed = SEED,
    quiet = TRUE
  )
  expect_s4_class(fit, "lgpfit")
  si <- get_stan_input(fit)
  expect_equal(si$is_likelihood_skipped, 1)
  p1 <- plot_warp(fit)
  p2 <- plot_effect_times(fit, verbose = FALSE)
  p3 <- plot_pred(fit, alpha = 1.0, draws = c(3, 4, 6), verbose = FALSE)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("Model with heterogeneous disease effect (f sampled)", {
  formula <- y ~ gp(age) + het(id) * gp_ns(dis_age)
  data <- testdata_001
  data$y <- round(exp(data$y))
  SW({
    fit <- lgp(
      formula = formula,
      data = data,
      prior_only = TRUE,
      prior = list(wrp = igam(14, 5)),
      likelihood = "nb",
      iter = N_ITER,
      chains = N_CHAINS,
      refresh = 0,
      seed = SEED,
      quiet = TRUE
    )
  })

  expect_s4_class(fit, "lgpfit")
  si <- get_stan_input(fit)
  expect_equal(si$is_likelihood_skipped, 1)
  p1 <- plot_beta(fit, verbose = FALSE)
  expect_s3_class(p1, "ggplot")
})
