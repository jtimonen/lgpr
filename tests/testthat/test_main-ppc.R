library(lgpr)
source("helpers/SW.R")

ITER <- 15
CHAINS <- 1
SEED <- 287

# -------------------------------------------------------------------------

context("ppc() function")


test_that("ppc() works for Gaussian model", {
  SW({
    fit <- example_fit(
      sample_f = TRUE,
      iter = ITER,
      chains = CHAINS,
      seed = SEED,
      prior_only = TRUE
    )
  })
  plt <- ppc(fit)
  expect_s3_class(plt, "ggplot")
})

test_that("ppc() works for Poisson model", {
  SW({
    fit <- example_fit(
      likelihood = "Poisson",
      iter = ITER,
      chains = CHAINS,
      seed = SEED
    )
  })
  plt <- ppc(fit)
  expect_s3_class(plt, "ggplot")
})

test_that("ppc() works for NB model", {
  SW({
    fit <- example_fit(
      likelihood = "NB",
      iter = ITER,
      chains = CHAINS,
      seed = SEED
    )
  })
  plt <- ppc(fit)
  expect_s3_class(plt, "ggplot")
})

test_that("ppc() works for binomial model", {
  SW({
    fit <- example_fit(
      likelihood = "binomial",
      iter = ITER,
      chains = CHAINS,
      seed = SEED,
      prior_only = TRUE
    )
  })
  plt <- ppc(fit)
  expect_s3_class(plt, "ggplot")
})


test_that("ppc() works for beta-binomial model", {
  SW({
    fit <- example_fit(
      likelihood = "bb",
      iter = ITER,
      chains = CHAINS,
      seed = SEED
    )
  })
  plt <- ppc(fit)
  expect_s3_class(plt, "ggplot")
})
