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
  a <- parse_formula(as.formula("y ~ gp(x) + zerosum(a)"))
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
model <- create_model(y ~ uncrt(id) * heter(id) * gp_warp(dis_age) +
  zerosum(sex) * gp(age),
prior = list(effect_time_info = et),
data = testdata_001
)

test_that("lgpmodel getters work", {
  n <- get_num_obs(model)
  df <- get_component_info(model)
  nams <- get_covariate_names(model)
  om <- get_obs_model(model)
  expect_equal(n, 24)
  expect_equal(dim(df), c(2, 9))
  expect_equal(nams, "dis_age, age, id, sex")
  expect_equal(om, "gaussian")
})

test_that("get_y works", {
  y1 <- get_y(model)
  y2 <- get_y(model, original = FALSE)
  expect_equal(length(y1), 24)
  expect_equal(length(y2), 24)
  expect_gt(stats::sd(y1), stats::sd(y2))
})

test_that("get_x_cat and get_x_cont work", {
  x <- get_x_cat(model)
  expect_equal(rownames(x), c("id", "sex"))
  x <- get_x_cont(model)
  expect_equal(sum(is.nan(x)), 12)
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

test_that("print_stan_input prints output", {
  expect_output(print_stan_input(model))
})

test_that("create_plot_df works correctly", {
  df <- create_plot_df(model)
  expect_equal(names(df), c("id", "age", "y"))
  a <- c(1,2,3)
  b <- seq(1, 24, by = 1)
  expect_error(create_plot_df(model, x = a), "length should be")
  reason <- "group_by has invalid type 'numeric'. Allowed types are"
  expect_error(create_plot_df(model, x = b, group_by = b), reason)
  ID <- as.factor(paste0("M", rep(c(1:4), each = 6)))
  df <- create_plot_df(model, x = b, group_by = ID)
  expect_equal(names(df), c("ID", "b", "y"))
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
  p <- plot_sim(dat, i_test = c(1, 2, 3), ncol = 4)
  expect_s3_class(p, "ggplot")
})

test_that("simulated data with disease effect can be plotted", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  p <- plot_sim(dat, i_test = c(1, 2, 3), ncol = 4)
  expect_s3_class(p, "ggplot")
})

test_that("show method for simualated data prints output", {
  dat <- simulate_data(
    N = 4,
    t_data = seq(6, 36, by = 6),
    covariates = c(0, 2),
    t_observed = "after_1"
  )
  expect_output(show(dat))
})


# -------------------------------------------------------------------------

context("Methods for lgpfit objects")

sim <- simulate_data(
  N = 4,
  t_data = seq(6, 36, by = 6),
  covariates = c(0, 2),
  t_observed = "after_1"
)
data <- sim@data
data$y <- data$y + 5

test_that("model fit can be visualized with data on original scale", {
  suppressWarnings({
    fit <- lgp(y ~ gp(age) + zerosum(id) + gp_warp_vm(diseaseAge),
      data = data, chains = 1, iter = 300, refresh = 0
    )
    p1 <- plot_fit(fit)
    p2 <- plot_fit(fit, draws = c(2, 3))
    p3 <- plot_fit(fit, draws = 3)
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    expect_s3_class(p3, "ggplot")
    rel <- relevances(fit)
    expect_equal(dim(rel), c(4, 2))
  })
})

test_that("fit summary prints output", {
  suppressWarnings({
    fit <- lgp(y ~ gp(age) + zerosum(z) * gp(age),
      data = data, chains = 1, iter = 100, refresh = 0
    )
    expect_output(fit_summary(fit))
    expect_output(show(fit))
  })
})
