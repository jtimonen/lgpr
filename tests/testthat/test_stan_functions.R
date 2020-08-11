library(lgpr)
library(rstan)

# stan_files is a symlink to <pkgdir>/inst/stan/
HOME <- file.path("stan_files", "chunks")
FILES <- c(
  "functions-utils.stan",
  "functions-kernels_base.stan",
  "functions-kernels_single.stan",
  "functions-kernels_many.stan",
  "functions-posterior.stan",
  "functions-prior.stan"
)

# Create stan model containing only the functions block
f_list <- lapply(file.path(HOME, FILES), FUN = readLines)
functions <- paste(unlist(f_list), collapse = "\n")
model_code <- paste(c("functions {", functions, "}"), collapse = "\n")

# Build the stan model
stanc_ret <- rstan::stanc(
  model_code = model_code,
  model_name = "only_functions",
  allow_undefined = TRUE
)

# Expose the functions
sf <- rstan::expose_stan_functions(stanc_ret,
  rebuild = TRUE,
  verbose = TRUE,
  show_compiler_warnings = FALSE
)

# 0. SETUP ----------------------------------------------------------------

context("Stan functions: setup")

test_that("Stan files are found", {
  expect_true(file.exists(HOME))
})

test_that("there is a correct number of exposed Stan functions", {
  expect_equal(length(sf), 21)
})


# 1. STAN UTILS -----------------------------------------------------------

context("Stan utils: input warping")

test_that("STAN_warp_input works for scalar input", {
  w <- STAN_warp_input(-1, 1.32)
  expect_equal(w, -0.5783634)
})

test_that("STAN_warp_input works for vector input", {
  w <- STAN_warp_input(c(-1, 0, 1), 1.32)
  w_correct <- c(-0.5783634, 0.0, 0.5783634)
  expect_equal(w, w_correct)
})

test_that("STAN_warp_input errors with invalid steepness input", {
  expect_error(STAN_warp_input(1, -1))
  expect_error(STAN_warp_input(1, 0.0))
  expect_error(STAN_warp_input(1, NaN))
  expect_error(STAN_warp_input(1, Inf))
})

test_that("STAN_warp_input works similarly as lgpr:::warp_input", {
  a <- exp(stats::rnorm(1)) # random steepness
  x <- seq(-3, 3, by = 1.33)
  expect_equal(
    STAN_warp_input(x, a),
    lgpr:::warp_input(x, a, 0, 1)
  )
})


context("Stan utils: variance masking")

test_that("STAN_var_mask works for scalar input", {
  m <- STAN_var_mask(-1, 1.32)
  expect_equal(m, 0.2108183)
})

test_that("STAN_var_mask works for vector input", {
  m <- STAN_var_mask(c(-1, 0, 1), 1.32)
  m_correct <- c(0.2108183, 0.5, 0.7891817)
  expect_equal(m, m_correct)
})

test_that("STAN_var_mask errors with invalid steepness input", {
  expect_error(STAN_var_mask(1, -1))
  expect_error(STAN_var_mask(1, 0.0))
  expect_error(STAN_var_mask(1, NaN))
  expect_error(STAN_var_mask(1, Inf))
})

test_that("STAN_var_mask works similarly as lgpr:::var_mask", {
  a <- 0.6
  x <- c(-5, 0, 5)
  expect_equal(
    STAN_var_mask(x, a),
    lgpr:::var_mask(x, a)
  )
})


context("Stan utils: expanding a vector")

test_that("STAN_expand works for valid input", {
  p <- c(0.1, 0.2)
  v <- STAN_expand(p, c(2, 3, 2, 3))
  expect_equal(v, c(0.1, 0.2, 0.1, 0.2))
})

test_that("STAN_expand errors when idx_expand has out of bounds indices", {
  p <- c(0.1, 0.2)
  idx1 <- c(2, 3, 0, 3)
  idx2 <- c(2, 3, 4, 3)
  expect_error(STAN_expand(p, idx1))
  expect_error(STAN_expand(p, idx2))
})


context("Stan utils: editing disease-related age")

test_that("STAN_edit_dis_age works properly", {
  x_dis_age <- c(
    -24, -12, 0, 12, -24, -12, 0, 12,
    0, 0, 0, 0, 0, -12, 0, 12, 16
  )
  teff <- c(-1, 2, 10)
  teff_obs <- c(0, 6, 12)
  case_ids <- c(1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 3, 3, 3, 3)
  idx_expand <- case_ids + 1
  expand_expect <- c(-1, -1, -1, -1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 10, 10, 10, 10)
  expect_equal(STAN_expand(teff, idx_expand), expand_expect)
  t_edit <- STAN_edit_dis_age(x_dis_age, idx_expand, teff_obs, teff)
  t_expect <- c(-23, -11, 1, 13, -20, -8, 4, 16, 0, 0, 0, 0, 0, -10, 2, 14, 18)
  expect_equal(t_edit, t_expect)
})

# 2. STAN BASE KERNELS ----------------------------------------------------

context("Stan base kernels: zero-sum kernel")

test_that("zero-sum kernel works correctly", {
  M <- 2
  x <- c(1, 1, 2)
  a <- c(
    1, 1, -1,
    1, 1, -1,
    -1, -1, 1
  )
  K_expect <- matrix(a, 3, 3, byrow = TRUE)
  K <- STAN_kernel_base_zerosum(x, x, M)
  expect_equal(K, K_expect)
})

test_that("zero-sum kernel works similarly in R and Stan", {
  M <- 3
  x <- sample.int(M, size = 8, replace = TRUE)
  expect_equal(
    STAN_kernel_base_zerosum(x, x, M),
    kernel_zerosum(x, x, M)
  )
})

test_that("zero-sum kernel errors if number of categories is one", {
  M <- 1
  x <- sample.int(M, size = 8, replace = TRUE)
  expect_error(STAN_kernel_base_zerosum(x, x, M))
})


context("Stan base kernels: ordinary categorical kernel")

test_that("categorical kernel works correctly", {
  x <- c(1, 1, 2)
  a <- c(
    1, 1, 0,
    1, 1, 0,
    0, 0, 1
  )
  K_expect <- matrix(a, 3, 3, byrow = TRUE)
  K <- STAN_kernel_base_cat(x, x)
  expect_equal(K, K_expect)
})


context("Stan base kernels: binary mask kernel")

test_that("binary mask kernel works correctly", {
  x <- c(1, 1, 2)
  a <- c(
    1, 1, 0,
    1, 1, 0,
    0, 0, 0
  )
  b <- c(
    0, 0, 0,
    0, 0, 0,
    0, 0, 1
  )
  K1_expect <- matrix(a, 3, 3, byrow = TRUE)
  K2_expect <- matrix(b, 3, 3, byrow = TRUE)
  K1 <- STAN_kernel_base_bin(x, x, 1)
  K2 <- STAN_kernel_base_bin(x, x, 2)
  expect_equal(K1, K1_expect)
  expect_equal(K2, K2_expect)
})

test_that("binary mask kernel works similarly in R and Stan", {
  x <- sample.int(2, size = 8, replace = TRUE) - 1
  expect_equal(
    STAN_kernel_base_bin(x, x, 1),
    kernel_bin(x, x)
  )
})


context("Stan base kernels: variance mask kernel")

test_that("variance mask kernel works correctly", {
  x <- c(12, 0, 12)
  stp <- 0.2
  vm_params <- c(0.05, 0.6)
  v <- c(
    0.9755191, 0.9382995, 0.9755191,
    0.9382995, 0.9025000, 0.9382995,
    0.9755191, 0.9382995, 0.9755191
  )
  K_expect <- matrix(v, 3, 3, byrow = TRUE)
  K <- STAN_kernel_base_var_mask(x, x, stp, vm_params)
  expect_equal(K, K_expect)
})


context("Stan base kernels: disease mask kernel")

test_that("disease mask kernel works correctly", {
  x1 <- c(1, 3, 0)
  x2 <- c(1, 0, 2, 0)
  a <- c(
    1, 0, 1, 0,
    1, 0, 1, 0,
    0, 0, 0, 0
  )
  K_expect <- matrix(a, 3, 4, byrow = TRUE)
  K <- STAN_kernel_base_disease_mask(x1, x2)
  expect_equal(K, K_expect)
})

test_that("variance mask kernel works similarly in R and Stan", {
  x <- c(-24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12)
  stp <- 1.0
  vm_params <- c(0.05, 0.6)
  K_stan <- STAN_kernel_base_var_mask(x, x, stp, vm_params)
  K_r <- kernel_var_mask(x, x, vm_params, stp)
  expect_equal(K_stan, K_r)
})

test_that("variance mask kernel errors if steepness is not valid", {
  x <- c(-24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12)
  vm_params <- c(0.05, 0.6)
  expect_error(STAN_kernel_base_var_mask(x, x, 0, vm_params))
  expect_error(STAN_kernel_base_var_mask(x, x, -1.0, vm_params))
  expect_error(STAN_kernel_base_var_mask(x, x, NaN, vm_params))
  expect_error(STAN_kernel_base_var_mask(x, x, Inf, vm_params))
})

test_that("variance mask kernel errors if <vm_params> is not valid", {
  x <- c(-24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12, -24, 12, 0, 12)
  vm_params <- c(0.05, 0.6)
  expect_error(STAN_kernel_base_var_mask(x, x, 1, c(-0.05, 0.6)))
  expect_error(STAN_kernel_base_var_mask(x, x, 1, c(1.6, 0.6)))
  expect_error(STAN_kernel_base_var_mask(x, x, 1, c(0.05, NaN)))
  expect_error(STAN_kernel_base_var_mask(x, x, 1, c(NaN, 0.6)))
})

# 3. STAN KERNEL ARRAYS ---------------------------------------------------

context("Stan kernels: fixed kernel arrays")

test_that("correct number of matrices is returned", {
  dat <- lgpr:::test_data_x(3)
  n1 <- length(dat$x1_disc[[1]])
  n2 <- length(dat$x2_disc[[1]])
  KF <- STAN_kernel_fixed_all(
    n1, n2, dat$x1_disc, dat$x2_disc,
    dat$num_levels, dat$components
  )
  expect_equal(length(KF), 6)
})

test_that("matrices of correct size are returned", {
  N <- 5
  dat <- lgpr:::test_data_x(N)
  n1 <- length(dat$x1_disc[[1]])
  n2 <- length(dat$x2_disc[[1]])
  KF <- STAN_kernel_fixed_all(
    n1, n2, dat$x1_disc, dat$x2_disc,
    dat$num_levels, dat$components
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

  KF <- STAN_kernel_fixed_all(
    4, 4, dat$x2_disc, dat$x2_disc,
    dat$num_levels, dat$components
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
  n1 <- length(dat$x1_disc[[1]])
  n2 <- length(dat$x2_disc[[1]])
  KF <- STAN_kernel_fixed_all(
    n1, n2, dat$x1_disc, dat$x2_disc,
    dat$num_levels, dat$components
  )
  alpha <- c(1, 1, 1, 1, 1, 1)
  ell <- c(1, 1, 1, 1, 1)
  x1 <- dat$x1_cont
  x2 <- dat$x2_cont
  K <- STAN_kernel_all(
    n1, n2, KF, dat$components, x1, x2,
    alpha, ell, 0.1, list(),
    list(), list(), list(), list(), list()
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
  KF <- STAN_kernel_fixed_all(
    n1, n1, dat$x1_disc, dat$x1_disc,
    dat$num_levels, dat$components
  )
  alpha <- 2 * c(1, 1, 1, 1, 1, 1)
  ell <- 12 * c(1, 1, 1, 1, 1)
  x1 <- dat$x1_cont
  K <- STAN_kernel_all(
    n1, n1, KF, dat$components, x1, x1,
    alpha, ell, 0.1, list(),
    list(), list(), list(), list(), list()
  )
  diff <- abs(K[[4]][2, 3] - 2.426123)
  expect_lt(diff, 1e-6)
})

# 4. STAN GP POSTERIOR ----------------------------------------------------

context("Stan GP posterior")

test_that("componentwise means sum to total mean", {
  N <- 3
  dat <- lgpr:::test_data_x(N)
  n1 <- length(dat$x1_disc[[1]])
  KF <- STAN_kernel_fixed_all(
    n1, n1, dat$x1_disc, dat$x1_disc,
    dat$num_levels, dat$components
  )
  alpha <- 2 * c(1, 1, 1, 1, 1, 1)
  ell <- 12 * c(1, 1, 1, 1, 1)
  x1 <- dat$x1_cont
  KX <- STAN_kernel_all(
    n1, n1, KF, dat$components, x1, x1,
    alpha, ell, 0.1, list(),
    list(), list(), list(), list(), list()
  )
  y <- rep(1, n1)

  fp <- STAN_gp_posterior(KX, y, 1e-6, 1.0)
  f_sum <- fp[[1]]
  for (j in 2:6) {
    f_sum <- f_sum + fp[[j]]
  }
  diff <- f_sum - fp[[7]]
  expect_lt(max(abs(diff)), 1e-6)
})

# 5. STAN PRIORS ----------------------------------------------------------

require(stats)

context("Stan priors: log density")

test_that("normal prior is correct", {
  x <- 0.333
  mu <- -0.11
  sigma <- 0.23
  log_p <- STAN_log_prior(x, c(2, 0), c(mu, sigma, 0))
  expect_equal(log_p, stats::dnorm(!!x, !!mu, !!sigma, log = TRUE))
})

test_that("student-t prior is correct", {
  x <- 0.333
  nu <- 16
  sigma <- 1
  log_p <- STAN_log_prior(x, c(3, 0), c(nu, sigma, 0))
  expect_equal(log_p, stats::dt(!!x, !!nu, log = TRUE))
})

test_that("gamma prior is correct", {
  x <- 0.333
  a <- 5
  b <- 8
  log_p <- STAN_log_prior(x, c(4, 0), c(a, b, 0))
  expect_equal(log_p, stats::dgamma(!!x, shape = !!a, rate = !!b, log = TRUE))
})

test_that("inverse-gamma prior is correct", {
  x <- 0.333
  a <- 5
  b <- 8
  log_p <- STAN_log_prior(x, c(5, 0), c(a, b, 0))
  expect_equal(
    log_p,
    lgpr:::dinvgamma_stanlike(!!x, alpha = !!a, beta = !!b, log = TRUE)
  )
})

test_that("log-normal prior is correct", {
  x <- 0.333
  mu <- -0.11
  sigma <- 0.23
  log_p <- STAN_log_prior(x, c(6, 0), c(mu, sigma, 0))
  expect_equal(log_p, stats::dlnorm(!!x, !!mu, !!sigma, log = TRUE))
})

context("Stan priors: log det of Jacobian of a transformation")

test_that("square transform is taken into account", {
  x <- 0.333
  mu <- -0.11
  sigma <- 0.23
  log_p <- STAN_log_prior(x, c(6, 1), c(mu, sigma, 0))
  expect_lt(
    log_p, # -38.9
    stats::dlnorm((!!x)^2, !!mu, !!sigma, log = TRUE) # -38.5
  )
})

context("Stan priors: error handling")

test_that("an error is thrown if <types> is misspecified", {
  x <- 1
  expect_error(STAN_log_prior(x, c(6), c(1, 1, 0)))
  expect_error(STAN_log_prior(x, c(0, 1), c(1, 1, 0)))
  expect_error(STAN_log_prior(x, c(7, 1), c(1, 1, 0)))
  expect_error(STAN_log_prior(x, c(1, 3), c(1, 1, 0)))
})

test_that("an error is thrown if <hyper> is misspecified", {
  x <- 1
  expect_error(STAN_log_prior(x, c(2, 0), c(1, 1)))
  expect_error(STAN_log_prior(x, c(2, 0), c(1, 1, 0, 1)))
  expect_error(STAN_log_prior(x, c(2, 0), c(1, -1, 0)))
})
