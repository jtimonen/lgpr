library(lgpr)

# -------------------------------------------------------------------------

context("Simulating data (sim)")

test_that("gaussian data can be simulated", {
  dat <- simulate_data(
    N = 8,
    t_data = seq(6, 36, by = 6),
    covariates = c(2, 2),
    lengthscales = c(6, 6, 6, 6),
    relevances = c(1, 1, 1, 0),
    names = c("sex", "location"),
    t_jitter = 0.5
  )
  data <- dat@data
  data_names <- c("id", "age", "sex", "location", "y")
  expect_equal(names(data), !!data_names)
  teff_names <- c("true", "observed")
  expect_equal(names(dat@effect_times), !!teff_names)
  info_names <- c("par_ell", "par_cont", "p_signal")
  expect_equal(names(dat@info), !!info_names)
})

test_that("poisson data can be simulated", {
  dat <- simulate_data(
    N = 8,
    t_data = seq(6, 36, by = 6),
    covariates = c(2, 2),
    noise_type = "poisson"
  )
  y <- dat@data$y
  diff <- round(y) - y
  expect_lt(max(diff), 1e-6)
})


test_that("negative binomial data can be simulated", {
  dat <- simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2, 2),
    noise_type = "nb"
  )
  y <- dat@data$y
  diff <- round(y) - y
  expect_lt(max(diff), 1e-6)
})

test_that("binomial data can be simulated", {
  dat <- simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2, 2),
    noise_type = "binomial",
    N_trials = 100
  )
  y <- dat@data$y
  diff <- round(y) - y
  expect_lt(max(diff), 1e-6)
})

test_that("different types of components can be simulated", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 1, 3, 4)
  )
  cnames <- c(
    "id.age", "age", "diseaseAge", "x", "offset", "group",
    "f", "g", "noise", "y"
  )
  expect_equal(names(dat@components), cnames)
})

test_that("disease component can be simulated only for selected individuals", {
  dat <- simulate_data(
    N = 8,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    N_affected = 2
  )
  n <- dim(dat@data)[1]
  s <- sum(abs(dat@components$diseaseAge[13:n]))
  expect_lt(s, 1e-6)
})


test_that("simulated effect time can be observed exactly", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "exact"
  )
  t <- c(21, 21, NaN, NaN)
  names(t) <- 1:4
  et <- dat@effect_times
  expect_equal(et$true, t)
  expect_equal(et$observed, t)
})


test_that("simulated effect time can be observed late", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  t1 <- c(21, 21, NaN, NaN)
  t2 <- c(30, 30, NaN, NaN)
  names(t1) <- 1:4
  names(t2) <- 1:4
  et <- dat@effect_times
  expect_equal(et$true, t1)
  expect_equal(et$observed, t2)
})

test_that("simulated effect time can be observed late randomly", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "random_0.5"
  )
  et <- dat@effect_times
  expect_equal(length(et$true), 4)
  expect_equal(length(et$observed), 4)
})

test_that("custom dis_fun argument can be given", {
  dis_fun <- function(x) {
    as.numeric(x > 0) + as.numeric(x < 20)
  }
  dat <- simulate_data(
    N = 4,
    t_effect_range = c(8, 10),
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    relevances = c(0, 0, 1, 0),
    dis_fun = dis_fun
  )
  ker <- dat@kernel_matrices
  expect_equal(dim(ker), c(24, 24, 4))
})


# -------------------------------------------------------------------------

context("Kernel functions (sim)")

test_that("base kernels work correctly", {
  expect_equal(
    sim_kernel_zerosum(c(1, 2), c(3, 2, 1), M = 3, alpha = 1),
    matrix(c(-0.5, -0.5, 1.0, -0.5, 1.0, -0.5),
      nrow = 2, ncol = 3, byrow = TRUE
    )
  )
  expect_equal(
    sim_kernel_bin(c(1, 2), c(3, 2, 1), pos_class = 2),
    matrix(c(0, 0, 0, 1, 0, 0), nrow = 2, ncol = 3, byrow = FALSE)
  )
  expect_equal(
    sim_kernel_se(-2, -2, ell = 20),
    matrix(1)
  )
  expect_equal(
    dim(sim_kernel_ns(c(1, 1, 2), c(0, 1), ell = 1, a = 1, b = -10, c = 1)),
    c(3, 2)
  )
})

test_that("sim_kernel_beta works correctly", {
  K <- sim_kernel_beta(c(0.1, 0.5, 1.0), c(1, 2, 3), c(1, 1, 2, 2, 3, 3))
  expect_equal(dim(K), c(3, 6))
  expect_equal(K[1, 1], 0.1)
  expect_equal(K[2, 4], 0.5)
  expect_equal(K[3, 6], 1.0)
})

test_that("base kernels give errors when supposed to", {
  expect_error(sim_kernel_se(0, c(-1, 0, 1), ell = 0))
  expect_error(sim_kernel_cat(0, c(-1, 0, 1), ell = -3, alpha = 1))
  expect_error(sim_kernel_ns(0, c(-1, 0, 1),
    ell = 1, alpha = -1, a = 1, b = 0, c = 1
  ))
  expect_error(sim_kernel_ns(0, c(-1, 0, 1),
    alpha = 1, a = 1, b = 0, c = 1
  ))
})
