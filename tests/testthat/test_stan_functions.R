library(lgpr)

# 1. STAN UTILS -----------------------------------------------------------

context("Stan utils")

STREAM <- get_stream()
test_that("rstan::get_stream is in namespace", {
  expect_true(class(STREAM) == "externalptr")
})

test_that("STAN_warp_input works for scalar input", {
  w <- warp_input(-1, 1.32)
  expect_equal(w, -0.5783634)
})

test_that("STAN_warp_input works for vector input", {
  w <- warp_input(c(-1, 0, 1), 1.32)
  w_correct <- c(-0.5783634, 0.0, 0.5783634)
  expect_equal(w, w_correct)
})

test_that("STAN_warp_input errors with invalid steepness input", {
  expect_error(warp_input(1, -1))
  expect_error(warp_input(1, 0.0))
  expect_error(warp_input(1, NaN))
  expect_error(warp_input(1, Inf))
})

test_that("STAN_warp_input works similarly as reference", {
  a <- exp(stats::rnorm(1)) # random steepness
  x <- seq(-3, 3, by = 1.33)
  expect_equal(
    warp_input(x, a),
    sim_warp_input(x, a, 0, 1)
  )
})

test_that("STAN_var_mask works for scalar input", {
  m <- var_mask(-1, 1.32)
  expect_equal(m, 0.2108183)
})

test_that("STAN_var_mask works for vector input", {
  m <- var_mask(c(-1, 0, 1), 1.32)
  m_correct <- c(0.2108183, 0.5, 0.7891817)
  expect_equal(m, m_correct)
})

test_that("STAN_var_mask errors with invalid steepness input", {
  expect_error(var_mask(1, -1))
  expect_error(var_mask(1, 0.0))
  expect_error(var_mask(1, NaN))
  expect_error(var_mask(1, Inf))
})

test_that("STAN_var_mask works similarly as reference", {
  a <- 0.6
  x <- c(-5, 0, 5)
  expect_equal(
    var_mask(x, a),
    sim_var_mask(x, a)
  )
})

test_that("STAN_expand works for valid input", {
  p <- c(0.1, 0.2)
  v <- STAN_expand(p, c(2, 3, 2, 3), STREAM)
  expect_equal(v, c(0.1, 0.2, 0.1, 0.2))
})

test_that("STAN_expand errors when idx_expand has out of bounds indices", {
  p <- c(0.1, 0.2)
  idx1 <- c(2, 3, 0, 3)
  idx2 <- c(2, 3, 4, 3)
  expect_error(STAN_expand(p, idx1, STREAM))
  expect_error(STAN_expand(p, idx2, STREAM))
})

test_that("STAN_edit_x_cont works properly", {
  x_dis_age <- c(
    -24, -12, 0, 12, -24, -12, 0, 12,
    0, 0, 0, 0, 0, -12, 0, 12, 16
  )
  teff <- c(-1, 2, 10)
  teff_obs <- c(0, 6, 12)
  case_ids <- c(1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 3, 3, 3, 3)
  idx_expand <- case_ids + 1
  expand_expect <- c(-1, -1, -1, -1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 10, 10, 10, 10)
  expect_equal(STAN_expand(teff, idx_expand, STREAM), expand_expect)
  t_edit <- STAN_edit_x_cont(x_dis_age, idx_expand, teff_obs, teff, STREAM)
  t_expect <- c(-23, -11, 1, 13, -20, -8, 4, 16, 0, 0, 0, 0, 0, -10, 2, 14, 18)
  expect_equal(t_edit, t_expect)
})

test_that("STAN_check_prob_positive works properly", {
  expect_error(STAN_check_prob_positive(1.1, STREAM))
  expect_error(STAN_check_prob_positive(-0.1, STREAM))
  expect_error(STAN_check_prob_positive(0, STREAM))
})

test_that("STAN_check_real_positive works properly", {
  expect_error(STAN_check_real_positive(-12, STREAM))
  expect_error(STAN_check_real_positive(0, STREAM))
})


# 2. STAN PRIORS ----------------------------------------------------------

require(stats)
context("Stan priors")

test_that("normal prior is correct", {
  x <- 0.333
  mu <- -0.11
  sigma <- 0.23
  log_p <- STAN_log_prior(x, c(2, 0), c(mu, sigma, 0), STREAM)
  expect_equal(log_p, stats::dnorm(!!x, !!mu, !!sigma, log = TRUE))
})

test_that("student-t prior is correct", {
  x <- 0.333
  nu <- 16
  sigma <- 1
  log_p <- STAN_log_prior(x, c(3, 0), c(nu, sigma, 0), STREAM)
  expect_equal(log_p, stats::dt(!!x, !!nu, log = TRUE))
})

test_that("gamma prior is correct", {
  x <- 0.333
  a <- 5
  b <- 8
  log_p <- STAN_log_prior(x, c(4, 0), c(a, b, 0), STREAM)
  expect_equal(log_p, stats::dgamma(!!x, shape = !!a, rate = !!b, log = TRUE))
})

test_that("inverse-gamma prior is correct", {
  x <- 0.333
  a <- 5
  b <- 8
  log_p <- STAN_log_prior(x, c(5, 0), c(a, b, 0), STREAM)
  expect_equal(
    log_p,
    lgpr:::dinvgamma_stanlike(!!x, alpha = !!a, beta = !!b, log = TRUE)
  )
})

test_that("log-normal prior is correct", {
  x <- 0.333
  mu <- -0.11
  sigma <- 0.23
  log_p <- STAN_log_prior(x, c(6, 0), c(mu, sigma, 0), STREAM)
  expect_equal(log_p, stats::dlnorm(!!x, !!mu, !!sigma, log = TRUE))
})

test_that("square transform is taken into account", {
  x <- 0.333
  mu <- -0.11
  sigma <- 0.23
  log_p <- STAN_log_prior(x, c(6, 1), c(mu, sigma, 0), STREAM)
  expect_lt(
    log_p, # -38.9
    stats::dlnorm((!!x)^2, !!mu, !!sigma, log = TRUE) # -38.5
  )
})
