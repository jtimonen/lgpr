library(lgpr)
source("helpers/SW.R")

ITER <- 17
CHAINS <- 1
SEED <- 939

# -------------------------------------------------------------------------

context("draw_pred() function")


test_that("draw_pred() works for Gaussian model", {
  SW({
    fit <- example_fit(
      sample_f = TRUE,
      iter = ITER,
      chains = CHAINS,
      seed = SEED
    )
  })
  y_draws <- draw_pred(fit)
  expect_equal(length(dim(y_draws)), 2)
  expect_equal(dim(y_draws)[2], 30)
})

test_that("draw_pred() works for Poisson model", {
  SW({
    fit <- example_fit(
      likelihood = "Poisson",
      iter = ITER,
      chains = CHAINS,
      seed = SEED
    )
  })
  y_draws <- draw_pred(fit)
  expect_equal(length(dim(y_draws)), 2)
  expect_equal(dim(y_draws)[2], 30)
})

test_that("draw_pred() works for NB model", {
  SW({
    fit <- example_fit(
      likelihood = "NB",
      iter = ITER,
      chains = CHAINS,
      seed = SEED
    )
  })
  y_draws <- draw_pred(fit)
  expect_equal(length(dim(y_draws)), 2)
  expect_equal(dim(y_draws)[2], 30)
})

test_that("draw_pred() works for binomial model", {
  SW({
    fit <- example_fit(
      likelihood = "binomial",
      iter = ITER,
      chains = CHAINS,
      seed = SEED
    )
  })
  y_draws <- draw_pred(fit)
  expect_equal(length(dim(y_draws)), 2)
  expect_equal(dim(y_draws)[2], 30)
})


test_that("draw_pred() works for beta-binomial model", {
  SW({
    fit <- example_fit(
      likelihood = "bb",
      iter = ITER,
      chains = CHAINS,
      seed = SEED
    )
  })
  y_draws <- draw_pred(fit)
  expect_equal(length(dim(y_draws)), 2)
  expect_equal(dim(y_draws)[2], 30)
})
