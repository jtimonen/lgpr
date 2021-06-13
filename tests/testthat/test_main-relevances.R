library(lgpr)
source("helpers/SW.R")

# -------------------------------------------------------------------------
N_ITER <- 28 # should be even
N_CHAINS <- 1

context("Computing relevances")

test_that("relevances can be computed when sample_f = FALSE", {
  SW({
    fit <- example_fit(iter = N_ITER, chains = N_CHAINS)
  })
  r_before <- relevances(fit)
  expect_equal(length(r_before), 5)
  expect_equal(sum(r_before), 1.0)

  # After clearing postproc info
  fit <- clear_postproc(fit)
  r_after <- relevances(fit, verbose = FALSE)
  expect_equal(r_before, r_after)
})

test_that("relevances can be computed when sample_f = TRUE", {
  SW({
    fit <- example_fit(iter = N_ITER, chains = N_CHAINS, likelihood = "Poisson")
  })
  r_before <- relevances(fit)
  expect_equal(length(r_before), 5)
  expect_equal(sum(r_before), 1.0)

  # After clearing postproc info
  fit <- clear_postproc(fit)
  r_after <- relevances(fit, verbose = FALSE)
  expect_equal(r_before, r_after)
})



test_that("relevances can be computed without reduce", {
  SW({
    fit <- example_fit(iter = N_ITER, chains = N_CHAINS)
  })
  r_before <- relevances(fit, reduce = NULL)
  expected_dim <- c(N_ITER / 2, 5)
  expect_equal(dim(r_before), expected_dim)
  is_one <- abs(rowSums(r_before) - 1) < 1e-6
  expect_equal(sum(!is_one), 0)

  # After clearing postproc info
  fit <- clear_postproc(fit)
  r_after <- relevances(fit, reduce = NULL, verbose = FALSE)
  expect_equal(dim(r_after), expected_dim)
  is_one <- abs(rowSums(r_after) - 1) < 1e-6
  expect_equal(sum(!is_one), 0)
})
