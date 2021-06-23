library(lgpr)

# 1. STAN UTILS -----------------------------------------------------------

context("STAN functions: utils")

STREAM <- get_stream()
test_that("rstan::get_stream is in namespace", {
  expect_true(class(STREAM) == "externalptr")
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

test_that("STAN_vectorsum works properly", {
  vecs <- list(c(1, 10, 1), c(2, 1, 4))
  s <- STAN_vectorsum(vecs, 3, STREAM)
  expect_equal(s, c(3, 11, 5))
})

# 2. STAN PRIORS ----------------------------------------------------------

require(stats)
context("STAN functions: priors")

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

n1 <- 323
test_that("input warping function works similarly in R and Stan code", {
  x <- sort(rnorm(n = n1))
  a <- 0.3 + runif(1)
  w1 <- warp_input(x, a)
  w2 <- STAN_warp_input(x, a, get_stream())
  expect_equal(w1, w2)
})
