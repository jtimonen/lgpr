library(lgpr)
source("helpers/SW.R")

ITER <- 121
CHAINS <- 3
SEED <- 232

# -------------------------------------------------------------------------

context("Prior sampling with sample_param_prior")

test_that("Model with uncertain effect time (f marginalized)", {
  formula <- y ~ gp(age) + unc(id) * gp_vm(dis_age)
  data <- testdata_001
  et <- list(backwards = FALSE, lower = 15, upper = 30, zero = 0)
  prior <- list(effect_time_info = et, wrp = igam(14, 5))
  model <- create_model(
    formula = formula,
    data = data,
    prior = prior
  )

  # NOTE: Not suppressing warnings here!
  fit <- sample_param_prior(model,
    iter = 2000,
    chains = 1,
    refresh = 0,
    seed = SEED,
    quiet = TRUE
  )
  expect_s4_class(fit, "stanfit")
})

test_that("Model with heterogeneous disease effect (f sampled)", {
  formula <- y ~ gp(age) + het(id) * gp_ns(dis_age)
  data <- testdata_001
  data$y <- round(exp(data$y))
  model <- create_model(
    formula = formula,
    data = data,
    prior = list(wrp = igam(14, 5)),
    likelihood = "Poisson"
  )

  # NOTE: Not suppressing warnings here!
  fit <- sample_param_prior(model,
    iter = 2000,
    chains = 1,
    refresh = 0,
    seed = SEED,
    quiet = TRUE
  )
  expect_s4_class(fit, "stanfit")
})

# -------------------------------------------------------------------------

context("Using prior_pred")


test_that("Model with uncertain effect time (f marginalized)", {
  formula <- y ~ gp(age) + unc(id) * gp_vm(dis_age)
  data <- testdata_001
  et <- list(backwards = FALSE, lower = 15, upper = 30, zero = 0)
  prior <- list(effect_time_info = et, wrp = igam(14, 5))
  model <- create_model(
    formula = formula,
    data = data,
    prior = prior
  )
  SW({
    pp <- prior_pred(model,
      iter = ITER, chains = CHAINS,
      seed = SEED, quiet = TRUE
    )
  })

  expect_equal(names(pp), c("y_draws", "pred_draws", "param_draws"))
})


test_that("Model with heterogeneous disease effect (f sampled)", {
  formula <- y ~ gp(age) + het(id) * gp_ns(dis_age)
  data <- testdata_001
  data$y <- round(exp(data$y))
  model <- create_model(
    formula = formula,
    data = data,
    prior = list(wrp = igam(14, 5)),
    likelihood = "NB"
  )
  SW({
    pp <- prior_pred(model, iter = ITER, chains = CHAINS, quiet = TRUE)
  })

  expect_equal(names(pp), c("y_draws", "pred_draws", "param_draws"))
})
