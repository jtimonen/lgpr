library(lgpr)
set.seed(123)

# -------------------------------------------------------------------------

context("Main function lgp")

# Create test input
sim <- simulate_data(
  N = 4,
  t_data = seq(6, 36, by = 6),
  covariates = c(0, 1, 2, 3),
  lengthscales = rep(12, 5),
  relevances = rep(1, 6),
  t_jitter = 0.5
)

# Formula
formula <- y ~ zerosum(id) * gp(age) + gp_warp_vm(diseaseAge) +
  categ(z) + gp(age) + gp(x)

test_that("an lgpfit object is returned and can be plotted", {
  suppressWarnings({
    fit <- lgp(
      formula = formula,
      data = sim@data,
      iter = 100,
      chains = 2,
      refresh = 0, cores = 2
    )
    expect_s4_class(fit, "lgpfit")
    p1a <- plot_posterior(fit)
    p1b <- plot_posterior(fit, type = "trace")
    p1c <- plot_posterior(fit, type = "dens", regex_pars = "alpha")
    p2 <- plot_posterior_warp(fit)
    expect_s3_class(p1a, "ggplot")
    expect_s3_class(p1b, "ggplot")
    expect_s3_class(p1c, "ggplot")
    expect_s3_class(p2, "ggplot")
  })
})

test_that("sbc can be performed", {
  suppressWarnings({
    fit <- calibrate(
      formula = formula,
      data = sim@data,
      iter = 1000,
      chains = 2,
      refresh = 0
    )
    expect_s4_class(fit, "lgpfit")
    p <- plot_posterior(fit)
    expect_s3_class(p, "ggplot")
  })
})


# Create test input
sim <- simulate_data(
  N = 4,
  t_data = seq(6, 36, by = 6),
  covariates = c(0),
  lengthscales = rep(12, 3),
  relevances = rep(1, 3),
  t_jitter = 0.5
)

# Formula
formula <- y ~ zerosum(id) * gp(age) + uncrt(id) * gp_warp_vm(diseaseAge)

test_that("a model with uncertain disease age needs prior specified", {
  data <- sim@data
  reason <- "you must specify 'effect_time_info' in the prior list"
  expect_error(create_model(formula = formula, data = data), reason)
})

test_that("a model with uncertain disease age can be created and fit", {
  et <- list(backwards = FALSE, lower = 15, upper = 30, zero = 0)
  prior <- list(effect_time_info = et)
  model <- create_model(formula = formula, data = sim@data, prior = prior)
  si <- get_stan_input(model)
  expect_equal(si$num_bt, 2)
  expect_equal(si$num_vm, 1)
  expect_equal(si$num_ns, 1)
  expect_equal(si$num_heter, 0)
  suppressWarnings({
    fit <- sample_model(model,
      iter = 400,
      chains = 2,
      refresh = 0,
      cores = 2
    )
    expect_s4_class(fit, "lgpfit")
  })
})
