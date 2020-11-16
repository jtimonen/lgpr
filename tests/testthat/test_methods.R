library(lgpr)

# -------------------------------------------------------------------------

context("Validation of expr, formula and scaling objects")

test_that("lgpexpr validation works correctly", {
  a <- lgpexpr(fun = "gp", covariate = "x")
  expect_true(check_lgpexpr(a))
  b <- a
  b@fun <- "moi"
  msg <- check_lgpexpr(b)
  expect_error(stop(msg), "<fun> must be one of")
  c <- a
  c@covariate <- ""
  msg <- check_lgpexpr(c)
  expect_error(stop(msg), "covariate name cannot be empty")
})

test_that("lgpformula validation works correctly", {
  a <- parse_formula(as.formula("y ~ gp(x) + zs(a)"), NULL)
  expect_true(check_lgpformula(a))
  b <- a
  b@y_name <- "x"
  msg <- check_lgpformula(b)
  expect_error(stop(msg), "response variable cannot be also")
})

test_that("lgpscaling validation works correctly", {
  f1 <- function(x) x / 2
  f2 <- function(x) 2 * x
  a <- new("lgpscaling", fun = f1, fun_inv = f2, var_name = "x")
  expect_true(check_lgpscaling(a))
  b <- a
  b@var_name <- ""
  msg <- check_lgpscaling(b)
  expect_error(stop(msg), "name length must be at least 1")
  c <- a
  c@fun <- function(x) 3 * x
  msg <- check_lgpscaling(c)
  expect_error(stop(msg), "<f_inv> is not an inverse function of <f>")
})

# -------------------------------------------------------------------------

context("Methods for lgpmodel objects")

et <- list(lower = 10, zero = 0, backwards = FALSE, upper = 20)
model <- create_model(y ~ unc(id) * het(id) * gp_ns(dis_age) +
  zs(sex) * gp(age),
prior = list(effect_time_info = et, wrp = igam(14, 5)),
data = testdata_001
)

test_that("lgpmodel getters work", {
  n <- get_num_obs(model)
  df <- component_info(model)
  om <- get_obs_model(model)
  expect_equal(n, 24)
  expect_equal(dim(df), c(2, 9))
  expect_equal(om, "gaussian")
})

test_that("get_y works", {
  y1 <- get_y(model)
  y2 <- get_y(model, original = FALSE)
  expect_equal(length(y1), 24)
  expect_equal(length(y2), 24)
  expect_gt(stats::sd(y1), stats::sd(y2))
})

test_that("getters work correctly for one-component models", {
  model <- create_model(y ~ age, data = testdata_001)
  expect_null(covariate_info.cat(model))
  model <- create_model(y ~ id, data = testdata_001)
  expect_null(covariate_info.cont(model))
})

test_that("prior summary works", {
  ps <- prior_summary(model, digits = 4)
  expect_equal(length(ps$Parameter), 9)
})

test_that("model summary prints output", {
  expect_output(model_summary(model))
  expect_output(show(model))
  expect_output(show(model@model_formula))
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
  expect_output(plot_sim(dat, verbose = TRUE))
  expect_output(plot_sim(dat, comp_idx = 1, verbose = TRUE))
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
