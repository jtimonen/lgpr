library(lgpr)
library(rstan)

context("Stan stream")
STREAM <- get_stream()

# -------------------------------------------------------------------------

context("Stan GP posterior")

test_that("componentwise means sum to total mean", {
  N <- 3
  dat <- lgpr:::test_data_x(N)
  n1 <- length(dat$x1_disc[[1]])
  KF <- STAN_kernel_const_all(
    n1, n1, dat$x1_cat, dat$x1_cat, dat$x1_cont_mask, dat$x1_cont_mask,
    dat$num_levels, dat$components, STREAM
  )
  alpha <- 2 * c(1, 1, 1, 1, 1, 1)
  ell <- 12 * c(1, 1, 1, 1, 1)
  x1 <- dat$x1_cont
  KX <- STAN_kernel_all(
    n1, n1, KF, dat$components, x1, x1,
    alpha, ell, 0.1, list(),
    list(), list(), list(), list(), list(), STREAM
  )
  y <- rep(1, n1)

  fp <- STAN_gp_posterior(KX, y, 1e-6, 1.0, STREAM)
  f_sum <- fp[[1]]
  for (j in 2:6) {
    f_sum <- f_sum + fp[[j]]
  }
  diff <- f_sum - fp[[7]]
  expect_lt(max(abs(diff)), 1e-6)
})
