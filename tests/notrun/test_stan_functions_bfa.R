library(lgpr)

context("Basis function approximation functions")
STREAM <- get_stream()

test_that("dirichlet eigenfunctions are computed correctly", {
  eig_fun <- function(x, m, L) {
    1 / sqrt(L) * sin((pi * m) / (2 * L) * (x + L))
  }
  x <- seq(-3, 3, by = 0.2)
  f1 <- STAN_phi(x, 1, 2.1, STREAM)
  f2 <- STAN_phi(x, 2, 1.1, STREAM)
  f3 <- STAN_phi(x, 3, 0.5, STREAM)
  expect_equal(sum(f1), sum(eig_fun(x, 1, 2.1)))
  expect_equal(sum(f2), sum(eig_fun(x, 2, 1.1)))
  expect_equal(sum(f3), sum(eig_fun(x, 3, 0.5)))
})

test_that("dirichlet eigenvalues are computed correctly", {
  eig_val <- function(m, L) {
    lam <- (pi * m) / (2 * L)
    return(lam^2)
  }
  l1 <- STAN_lambda(1, 2.1, STREAM)
  l2 <- STAN_lambda(2, 1.1, STREAM)
  expect_equal(l1, eig_val(1, 2.1))
  expect_equal(l2, eig_val(2, 1.1))
})

test_that("spectral density of exp quad kernel is correct", {
  spd <- function(w, l) {
    l * sqrt(2 * pi) * exp(-0.5 * (l * w)^2)
  }
  s1 <- STAN_spd_eq(0, 2, STREAM)
  s2 <- STAN_spd_eq(0.2, 1.5, STREAM)
  expect_equal(s1, spd(0, 2))
  expect_equal(s2, spd(0.2, 1.5))
})

test_that("multivariate normal approximation works", {
  n <- 3
  RM <- 5
  y <- rnorm(n)
  D_diag <- rep(2.3, RM)
  V <- matrix(rnorm(n * RM), n, RM)
  sigma <- 0.3
  p <- STAN_multi_normal_bfa_logpdf(y, V, D_diag, sigma, STREAM)
  expect_false(is.nan(p))
})
