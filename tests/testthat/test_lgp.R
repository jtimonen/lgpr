library(lgpr)

# Warnings are suppressed because sampling won't converge in such a short time

test_that("lgp runs", {
  expect_identical(
    suppressWarnings({
      as.vector(lgp(
        formula = y ~ id + age,
        data = simulate_data(N = 4, 10 * c(1, 2, 3, 4, 5))$data,
        iter = 100,
        chains = 1,
        refresh = 0,
        verbose = FALSE,
        relevance_method = "alpha"
      )@model@stan_dat$D)
    }),
    c(1, 1, 0, 0, 0, 0)
  )

  expect_identical(
    suppressWarnings({
      as.vector(lgp(
        formula = y ~ id + age + z,
        data = simulate_data(
          N = 4,
          t_data = 10 * c(1, 2, 3, 4, 5),
          covariates = c(2)
        )$data,
        iter = 100,
        chains = 1,
        verbose = FALSE,
        refresh = 0
      )@model@stan_dat$D)
    }),
    c(1, 1, 0, 0, 1, 0)
  )
})

test_that("lgp runs without id*age component", {
  expect_identical(
    suppressWarnings({
      as.vector(lgp(
        formula = y ~ age + z,
        data = simulate_data(
          N = 4,
          t_data = 10 * c(1, 2, 3, 4, 5),
          covariates = c(2)
        )$data,
        iter = 100,
        chains = 1,
        time_variable = "age",
        id_variable = "id",
        verbose = FALSE,
        refresh = 0
      )@model@stan_dat$D)
    }),
    c(0, 1, 0, 0, 1, 0)
  )
})

test_that("lgp runs without age*id component but with shared age", {
  expect_identical(
    suppressWarnings({
      as.vector(lgp(
        formula = y ~ age + z,
        data = simulate_data(
          N = 4,
          t_data = 10 * c(1, 2, 3, 4, 5),
          covariates = c(2)
        )$data,
        iter = 100,
        chains = 1,
        verbose = FALSE,
        refresh = 0
      )@model@stan_dat$D)
    }),
    c(0, 1, 0, 0, 1, 0)
  )
})


test_that("lgp can sample F", {
  expect_identical(
    suppressWarnings({
      as.vector(lgp(
        formula = y ~ id + age + diseaseAge + x + z + offset,
        data = simulate_data(
          N = 4,
          t_data = 10 * c(1, 2, 3, 4, 5),
          covariates = c(0, 1, 2, 3)
        )$data,
        iter = 100,
        chains = 1,
        refresh = 0,
        verbose = FALSE,
        offset_vars = c("offset"),
        sample_F = T
      )@model@stan_dat$D)
    }),
    c(1, 1, 1, 1, 1, 1)
  )
})

test_that(paste0(
  "lgp can be used without the vm kernel",
  " and with special disease modeling features"
), {
  expect_identical(
    suppressWarnings({
      as.vector(lgp(
        formula = y ~ id + age + diseaseAge + group,
        data = simulate_data(
          N = 4,
          t_data = 10 * c(1, 2, 3, 4, 5),
          covariates = c(0, 1, 2, 4),
          dis_fun = "gp_ns"
        )$data,
        iter = 100,
        chains = 1,
        refresh = 0,
        verbose = FALSE,
        offset_vars = c("group"),
        equal_effect = FALSE,
        uncertain_effect_time = TRUE,
        variance_mask = FALSE
      )@model@stan_dat$D)
    }),
    c(1, 1, 1, 0, 0, 1)
  )
})


test_that("lgp can sample from prior", {
  expect_equal(
    suppressWarnings({
      dim(lgp(
        formula = y ~ id + age + z, data = simulate_data(
          N = 4,
          t_data = 10 * c(1, 2, 3, 4, 5), covariates = c(2)
        )$data,
        iter = 100,
        chains = 1,
        refresh = 0, skip_postproc = TRUE, likelihood = "none"
      )@diagnostics)
    }),
    c(9, 3)
  )
})

test_that("lgp can be run using Poisson observation model", {
  expect_equal(suppressWarnings({
    dim(lgp(
      formula = y ~ id + age,
      data = simulate_data(
        N = 4, t_data = 10 * c(1, 2, 3, 4, 5),
        noise_type = "Poisson"
      )$data,
      likelihood = "Poisson",
      iter = 60,
      chains = 1,
      verbose = FALSE,
      refresh = 0
    )@relevances$samples)
  }), c(30, 3))
})


test_that("lgp can be run using NB observation model", {
  expect_equal(suppressWarnings({
    dim(lgp(
      formula = y ~ id + age,
      data = simulate_data(
        N = 4, t_data = 10 * c(1, 2, 3, 4, 5),
        noise_type = "NB",
        phi = 2
      )$data,
      likelihood = "NB",
      iter = 60,
      chains = 1,
      verbose = FALSE,
      refresh = 0
    )@relevances$samples)
  }), c(30, 3))
})


test_that("lgp can be run using binomial observation model", {
  expect_equal(suppressWarnings({
    dim(lgp(
      formula = y ~ id + age,
      data = simulate_data(
        N = 4, t_data = 10 * c(1, 2, 3, 4, 5),
        noise_type = "binomial", N_trials = 20
      )$data,
      likelihood = "binomial",
      N_trials = 20,
      iter = 60,
      chains = 1,
      verbose = FALSE,
      refresh = 0
    )@relevances$samples)
  }), c(30, 3))
})
