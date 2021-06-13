library(lgpr)

# -------------------------------------------------------------------------

context("Simulating data")

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
  info_names <- c("par_ell", "par_cont", "p_signal", "msg", "noise_type")
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
  exp_msg_contains <- "Line is the true signal mapped"
  expect_message(plot_sim(dat, verbose = TRUE), exp_msg_contains)
})

test_that("binomial data can be simulated", {
  dat <- simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2, 2),
    bin_kernel = TRUE,
    noise_type = "binomial",
    N_trials = 100
  )
  y <- dat@data$y
  diff <- round(y) - y
  expect_lt(max(diff), 1e-6)
})

test_that("beta-binomial data can be simulated", {
  dat <- simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2, 2),
    bin_kernel = TRUE,
    noise_type = "bb",
    N_trials = 100
  )
  y <- dat@data$y
  diff <- round(y) - y
  expect_lt(max(diff), 1e-6)
})

test_that("error is thrown with invalid input", {
  expect_error(simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2, 2),
    noise_type = "blaablaa",
  ), "given value 'blaablaa' for argument <likelihood> is invalid")

  expect_error(simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2),
    names = c("sex", "country")
  ), "lengths of names and covariates must match! found = 2 and 1")

  expect_error(simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2),
    relevances = c(1, 1, 1, 1)
  ), "must be")

  expect_error(simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2),
    n_categs = c(2, 3, 3)
  ), "The argument n_cat has invalid length")

  expect_error(simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2),
    n_categs = c(1)
  ), "must only contain integers larger than 1")

  expect_error(simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2),
    lengthscales = c(10, 10, 3, 3)
  ), "lengthscales has length 4, should be 3")

  expect_error(simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2),
    lengthscales = c(10, 10, 3, 3),
    N_affected = 1000
  ), "N_affected cannot be greater than")

  expect_error(simulate_data(
    N = 4,
    t_data = c(1, 2, 3),
    covariates = c(2, 1, 0),
    lengthscales = c(10, 10, 3, 3),
  ), "covariates vector must be increasing")
})

test_that("different types of components can be simulated", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 1, 3, 4),
    dis_fun = "gp_warp",
  )
  cnames <- c(
    "id.age", "age", "diseaseAge", "x", "offset", "group",
    "f", "h", "noise", "y"
  )
  expect_equal(names(dat@components), cnames)
})

test_that("invalid disease function cannot be given", {
  reason <- "dis_fun must be gp_warp or gp_warp_vm"
  expect_error(simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 1, 3, 4),
    dis_fun = "notavalidname",
  ), reason)
})

test_that("disease component can be simulated only for selected individuals", {
  dat <- simulate_data(
    N = 8,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = function(x) x + 0.1,
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


test_that("sim.check_too_far works correctly", {
  rem <- c(5, 7, 4)
  expect_true(sim.check_too_far(1, rem))
  reason <- "Not enough data points to go that far"
  expect_error(sim.check_too_far(4, rem), reason)
})



# -------------------------------------------------------------------------

context("Methods for lgpsim objects")
set.seed(123)

test_that("simulated data can be plotted", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(1, 2)
  )
  p1 <- plot_sim(dat, i_test = c(1, 2, 3), ncol = 4, verbose = FALSE)
  p2 <- plot_sim(dat, comp_idx = 2, verbose = FALSE) # not colored
  p3 <- plot_sim(dat, comp_idx = 3, color_by = "z", verbose = FALSE) # colored
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("simulated data with disease effect can be plotted", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  p1 <- plot_sim(dat, i_test = c(1, 2, 3), ncol = 4, verbose = FALSE) # vlines
  p2 <- plot_sim(dat, comp_idx = 1, verbose = FALSE)
  p3 <- plot_sim(dat, comp_idx = 3, color_by = "diseaseAge", verbose = FALSE)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_message(
    plot_sim(dat, verbose = TRUE),
    "Dashed vert. line is the 'observed' disease initiation time"
  )
  expect_message(
    plot_sim(dat, comp_idx = 2, verbose = TRUE),
    "component_idx = 2"
  )
})

test_that("show method for simulated data prints output", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  expect_output(show(dat))
})
