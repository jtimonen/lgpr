library(lgpr)
source("helpers/SW.R")
# source("tests/testthat/helpers/SW.R")

# -------------------------------------------------------------------------

context("Creating GaussianPrediction objects")

N_ITER <- 33
N_CHAINS <- 2
SW({
  fit <- example_fit(iter = N_ITER, chains = N_CHAINS)
})

test_that("pred works with defaults", {
  a <- pred(fit = fit, verbose = FALSE)
  expect_s4_class(a, "GaussianPrediction")
  expect_output(show(a))
})

test_that("pred works with reduce = NULL", {
  a <- pred(fit = fit, reduce = NULL, verbose = FALSE)
  expect_s4_class(a, "GaussianPrediction")
  expect_output(show(a))
})

test_that("pred works with defaults and given x", {
  x_pred <- new_x(fit, x_values = seq(0, 100, by = 2.5))
  a <- pred(fit = fit, x = x_pred, verbose = FALSE)
  expect_s4_class(a, "GaussianPrediction")
  expect_output(show(a))
})

test_that("pred works with defaults and given x, and reduce = NULL", {
  x_pred <- new_x(fit, x_values = seq(0, 100, by = 2.5))
  a <- pred(
    fit = fit, x = x_pred, reduce = NULL, draws = c(2, 5),
    verbose = FALSE
  )
  expect_s4_class(a, "GaussianPrediction")
  expect_output(show(a))

  # Check that default is verbose (should print also progbars)
  expect_output({
    pp <- pred(fit = fit, x = x_pred, reduce = NULL)
  })
})


# -------------------------------------------------------------------------

context("Creating Prediction objects")

SW({
  fit <- example_fit(iter = N_ITER, chains = N_CHAINS, likelihood = "nb")
})

test_that("pred works with defaults", {
  a <- pred(fit = fit, verbose = FALSE)
  expect_s4_class(a, "Prediction")
  expect_output(show(a))
})

test_that("pred works with reduce = NULL", {
  a <- pred(fit = fit, reduce = NULL, verbose = FALSE)
  expect_s4_class(a, "Prediction")
  expect_output(show(a))
})

test_that("pred works with defaults and given x", {
  x_pred <- new_x(fit, x_values = seq(0, 100, by = 2.5))
  a <- pred(fit = fit, x = x_pred, verbose = FALSE)
  expect_s4_class(a, "Prediction")
  expect_output(show(a))
})

test_that("pred works with defaults and given x, and reduce = NULL", {
  x_pred <- new_x(fit, x_values = seq(0, 100, by = 2.5))
  a <- pred(
    fit = fit, x = x_pred, reduce = NULL, draws = c(2, 5),
    verbose = FALSE
  )
  expect_s4_class(a, "Prediction")
  expect_output(show(a))
  
  # Check that default is verbose (should print also progbars)
  expect_output({
    pp <- pred(fit = fit, x = x_pred, reduce = NULL)
  })
})

