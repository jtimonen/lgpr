library(lgpr)
library(rstan)

context("Stan stream")
STREAM <- get_stream()

# -------------------------------------------------------------------------

context("Stan kernels: const kernel arrays")

test_that("correct number of matrices is returned", {
  dat <- lgpr:::test_data_x(3)
  n1 <- length(dat$x1_cat[[1]])
  n2 <- length(dat$x2_cat[[1]])
  KF <- STAN_kernel_const_all(
    n1, n2, dat$x1_cat, dat$x2_cat, dat$x1_cont_mask, dat$x2_cont_mask,
    dat$x_cat_num_levels, dat$components, STREAM
  )
  expect_equal(length(KF), 6)
})

test_that("matrices of correct size are returned", {
  N <- 5
  dat <- lgpr:::test_data_x(N)
  n1 <- length(dat$x1_disc[[1]])
  n2 <- length(dat$x2_disc[[1]])
  KF <- STAN_kernel_const_all(
    n1, n2, dat$x1_disc, dat$x2_disc,
    dat$num_levels, dat$components, STREAM
  )
  J <- length(KF)
  for (j in seq_len(J)) {
    expect_equal(dim(KF[[!!j]]), c(3 * N, 4))
  }
})

test_that("matrices with correct values are returned", {
  dat <- lgpr:::test_data_x(3)

  a1 <- c(
    1.0, 1.0, -0.5, -0.5,
    1.0, 1.0, -0.5, -0.5,
    -0.5, -0.5, 1.0, -0.5,
    -0.5, -0.5, -0.5, 1.0
  )

  a2 <- c(
    1.0, 1.0, -1.0, -1.0,
    1.0, 1.0, -1.0, -1.0,
    -1.0, -1.0, 1.0, 1.0,
    -1.0, -1.0, 1.0, 1.0
  )

  a3 <- c(
    1.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 1.0, 1.0
  )

  a4 <- rep(0, times = 16)
  a5 <- a4
  a5[11] <- 1.0
  a6 <- c(
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 1.0, 1.0
  )

  A1 <- matrix(a1, 4, 4, byrow = TRUE)
  A2 <- matrix(a2, 4, 4, byrow = TRUE)
  A3 <- matrix(a3, 4, 4, byrow = TRUE)
  A4 <- matrix(a4, 4, 4, byrow = TRUE)
  A5 <- matrix(a5, 4, 4, byrow = TRUE)
  A6 <- matrix(a6, 4, 4, byrow = TRUE)

  KF <- STAN_kernel_const_all(
    4, 4, dat$x2_cat, dat$x2_cat,
    dat$num_levels, dat$components, STREAM
  )
  expect_equal(KF[[1]], A1)
  expect_equal(KF[[2]], A2)
  expect_equal(KF[[3]], A3)
  expect_equal(KF[[4]], A4)
  expect_equal(KF[[5]], A5)
  expect_equal(KF[[6]], A6)
})

context("Stan kernels: full kernel array")

test_that("STAN_kernel_all can be used", {
  N <- 5
  dat <- lgpr:::test_data_x(N)
  n1 <- length(dat$x1_cat[[1]])
  n2 <- length(dat$x2_cat[[1]])
  KF <- STAN_kernel_const_all(
    n1, n2, dat$x1_cat, dat$x2_cat, dat$x1_cont_mask, dat$x2_cont_mask,
    dat$num_levels, dat$components, STREAM
  )
  alpha <- c(1, 1, 1, 1, 1, 1)
  ell <- c(1, 1, 1, 1, 1)
  x1 <- dat$x1_cont
  x2 <- dat$x2_cont
  K <- STAN_kernel_all(
    n1, n2, KF, dat$components, x1, x2,
    alpha, ell, 0.1, list(),
    list(), list(), list(), list(), list(), STREAM
  )
  J <- 6
  expect_equal(length(K), J)
  n1 <- length(x1[[1]])
  n2 <- length(x2[[1]])
  for (j in seq_len(J)) {
    expect_equal(dim(K[[!!j]]), c(!!n1, !!n2))
  }
})

test_that("STAN_kernel_all uses cov_exp_quad correctly", {
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
  K <- STAN_kernel_all(
    n1, n1, KF, dat$components, x1, x1,
    alpha, ell, 0.1, list(),
    list(), list(), list(), list(), list(), STREAM
  )
  diff <- abs(K[[4]][2, 3] - 2.426123)
  expect_lt(diff, 1e-6)
})
