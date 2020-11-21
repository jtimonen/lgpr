library(lgpr)

# -------------------------------------------------------------------------

context("Methods for lgpfit objects")

sim <- simulate_data(
  N = 4,
  t_data = seq(6, 36, by = 6),
  covariates = c(0, 2),
  t_observed = "after_1"
)

# Create first test fit object
data1 <- sim@data
data1$y <- data1$y + 5
suppressWarnings({
  fit1 <- lgp(y ~ gp(age) + zs(id) + gp_vm(diseaseAge),
    data = data1, chains = 1, iter = 100, refresh = 0
  )
})

# Create second test fit object
data2 <- data1
data2$y <- round(exp(data2$y - 5))
suppressWarnings({
  fit2 <- lgp(y ~ gp(age) + zs(id) + gp_vm(diseaseAge),
    likelihood = "nb",
    data = data2, chains = 1, iter = 100, refresh = 0
  )
  fit2_rng <- lgp(y ~ gp(age) + zs(id) + gp_vm(diseaseAge),
    likelihood = "nb",
    data = data2, chains = 1, iter = 100, refresh = 0,
    options = list(do_yrng = TRUE)
  )
})

test_that("fit summary prints output", {
  expect_output(fit_summary(fit1))
  expect_output(fit_summary(fit2))
  expect_output(show(fit1))
})

test_that("get_draws works correctly", {
  d <- get_draws(fit1)
  expect_equal(dim(d), c(50, 200))
  d <- get_draws(fit2)
  expect_equal(dim(d), c(50, 152))
  d <- get_draws(fit2_rng)
  expect_equal(dim(d), c(50, 152 + 24))
  d <- get_draws(fit1, pars = "f_post", include = FALSE)
  expect_equal(dim(d), c(50, 8))
  d <- get_draws(fit1, pars = "f_post", include = FALSE, reduce = mean)
  expect_equal(dim(d), c(1, 8))
  expect_equal(colnames(d)[1], "alpha[1]")
  d <- get_draws(fit1, pars = "f_post", include = FALSE, draws = 4)
  expect_equal(dim(d), c(1, 8))
  d <- get_draws(fit1, pars = "f_post", include = FALSE, draws = 11:20)
  expect_equal(dim(d), c(10, 8))
})


test_that("get_pred errors correctly", {
  reason <- "must be NULL"
  expect_error(get_pred(fit1, draws = 1, reduce = "mean"), reason)
})


test_that("get_pred.gaussian works correctly", {
  cnames <- c("gp(age)", "zs(id)", "gp_vm(diseaseAge)")
  a <- get_pred(fit1)
  expect_s4_class(a, "GaussianPrediction")
  expect_equal(names(a@f_comp_mean), cnames)
  expect_equal(names(a@f_comp_std), cnames)
  expect_equal(dim(a@y_std), c(50, 24))
  a <- get_pred(fit1, reduce = mean)
  expect_equal(dim(a@y_mean), c(1, 24))
  a <- get_pred(fit1, draws = c(14, 34))
  expect_equal(dim(a@f_std), c(2, 24))
  expect_true(check_GaussianPrediction(a))
  expect_output(show(a))
})

test_that("get_pred.sampled works correctly", {
  cnames <- c("gp(age)", "zs(id)", "gp_vm(diseaseAge)")
  a <- get_pred(fit2)
  expect_s4_class(a, "Prediction")
  expect_equal(names(a@f_comp), cnames)
  expect_equal(dim(a@h), c(50, 24))
  a <- get_pred(fit2, reduce = mean)
  expect_equal(dim(a@f), c(1, 24))
  a <- get_pred(fit2, draws = c(14, 34))
  expect_equal(dim(a@h), c(2, 24))
  expect_true(check_Prediction(a))
  expect_output(show(a))
})

test_that("relevances can be computed and sum to 1", {
  r1 <- relevances(fit1)
  r2 <- relevances(fit2)
  expect_equal(length(r1), 4)
  expect_equal(length(r2), 4)
  expect_equal(sum(r1), 1.0)
  expect_equal(sum(r2), 1.0)
})

test_that("select and select_freq work", {
  s1a <- select(fit1, reduce = mean)
  s1b <- select_freq(fit1)
  expect_equal(dim(s1a), c(4, 1))
  expect_equal(dim(s1b), c(4, 1))
  s2 <- select(fit2, threshold = 0.8)
  expect_equal(dim(s2), c(4, 1))
})

test_that("probabilistic selection works", {
  b <- select.integrate(fit1, reduce = stats::median, verbose = FALSE)
  expect_equal(dim(b$selected), c(101, 4))
  b <- select_freq.integrate(fit2, verbose = FALSE, h = 0.1)
  expect_equal(dim(b$freq), c(11, 4))
  expect_output(select.integrate(fit2, verbose = TRUE, h = 0.1))
  expect_output(select_freq.integrate(fit1, verbose = TRUE, h = 0.1))
})

test_that("predictions can be visualized with data on original scale", {
  p1 <- plot_pred(fit1, reduce = NULL)
  p2 <- plot_pred(fit1, reduce = NULL, draws = c(2, 3))
  p3 <- plot_pred(fit1, reduce = NULL, draws = 3)
  p4 <- plot_pred(fit1, reduce = mean)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  p1 <- plot_pred(fit2, reduce = NULL)
  p2 <- plot_pred(fit2, reduce = NULL, draws = c(2, 3))
  p3 <- plot_pred(fit2, reduce = NULL, draws = 3)
  p4 <- plot_pred(fit2, reduce = mean)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("inferred components can be visualized", {
  p1 <- plot_f(fit1, reduce = NULL)
  p2 <- plot_f(fit1,
    comp_idx = 3, reduce = NULL,
    draws = c(3, 4), color_by = "diseaseAge"
  )
  p3 <- plot_f(fit1,
    reduce = stats::median, color_by = "diseaseAge",
    alpha_err = 0.3, comp_idx = 2, MULT_STD = 0.03
  )
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  a <- plot_components(fit2,
    color_by = "diseaseAge", draw = FALSE,
    ylim = c(-3, 3), ncol = 2, nrow = 2
  )
  expect_s3_class(a[[1]], "ggplot")
  expect_equal(length(a), 4)
})

# -------------------------------------------------------------------------

context("Out-of-sample prediction")

t <- seq(0, 50, by = 5)
x1_pred <- new_x(data1, t, x_ns = "diseaseAge")
x2_pred <- new_x(data2, t, x_ns = "diseaseAge")

test_that("prediction kernel computations work correctly", {
  kers <- pred_kernels(fit1, x1_pred, NULL, NULL, FALSE)
  expect_equal(names(kers), c("data_vs_data", "pred_vs_data", "pred_vs_pred"))
  expect_equal(dim(kers$data_vs_data), c(50, 3, 24, 24))
  expect_equal(dim(kers$pred_vs_data), c(50, 3, 44, 24))
  expect_equal(dim(kers$pred_vs_pred), c(50, 3, 44, 44))
})

test_that("pred.gaussian works correctly", {
  out <- pred(fit1, x1_pred, verbose = FALSE)
  expect_s4_class(out, "GaussianPrediction")
  out <- pred(fit1, x1_pred, draws = c(3:43), verbose = FALSE)
  expect_s4_class(out, "GaussianPrediction")
  out <- pred(fit1, x1_pred, reduce = mean, verbose = FALSE)
  expect_s4_class(out, "GaussianPrediction")
  expect_equal(dim(out@y_mean), c(1, 44))
  expect_output(pred(fit1, x1_pred, verbose = TRUE))
  expect_output(pred(fit2, x2_pred, verbose = TRUE))
})

test_that("pred.kr works correctly", {
  out <- pred(fit2, x2_pred, reduce = mean, verbose = FALSE)
  expect_equal(dim(out@h), c(1, 44))
  out <- pred(fit2, x2_pred, reduce = NULL, verbose = FALSE)
  expect_equal(dim(out@h), c(50, 44))
})


# -------------------------------------------------------------------------

context("Out-of-sample prediction visualization")

test_that("out-of-sample predictions can be visualized", {
  os1 <- pred(fit1, x1_pred, verbose = FALSE)
  os2 <- pred(fit1, x1_pred,
    reduce = NULL, draws = c(8:10),
    verbose = FALSE
  )
  p1 <- plot_pred(fit1, pred = os1, x = x1_pred)
  p2 <- plot_pred(fit1, pred = os2, x = x1_pred)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("out-of-sample inferred components can be visualized", {
  os1 <- pred(fit1, x1_pred, verbose = FALSE)
  os2 <- pred(fit1, x1_pred,
    reduce = NULL, draws = c(8:10),
    verbose = FALSE
  )
  p1 <- plot_f(fit1, pred = os1, x = x1_pred)
  p2 <- plot_f(fit1,
    pred = os2, x = x1_pred, comp_idx = 3,
    color_by = "diseaseAge"
  )
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("out-of-sample KR can be visualized", {
  os1 <- pred(fit2, x2_pred, verbose = FALSE)
  os2 <- pred(fit2, x2_pred,
    reduce = NULL, draws = c(8:10),
    verbose = FALSE
  )
  p1 <- plot_pred(fit2, pred = os1, x = x2_pred)
  p2 <- plot_pred(fit2, pred = os2, x = x2_pred)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("out-of-sample KR components can be visualized", {
  os1 <- pred(fit2, x2_pred, verbose = FALSE)
  os2 <- pred(fit2, x2_pred,
    reduce = NULL, draws = c(8:10),
    verbose = FALSE
  )
  p1 <- plot_f(fit2, pred = os2, x = x2_pred)
  p2 <- plot_f(fit2,
    pred = os2, x = x2_pred, comp_idx = 3,
    color_by = "diseaseAge"
  )
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})
