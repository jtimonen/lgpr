library(lgpr)
source("helpers/SW.R")
# source("tests/testthat/helpers/SW.R")

# -------------------------------------------------------------------------
set.seed(4939)

context("Plotting GaussianPrediction objects")

N_ITER <- 42
N_CHAINS <- 1
SW({
  fit <- example_fit(iter = N_ITER, chains = N_CHAINS)
})

test_that("plotting works with defaults", {
  a <- pred(fit = fit, verbose = FALSE)
  p1 <- plot_pred(fit, pred = a)
  p2 <- plot_components(fit, pred = a, draw = FALSE)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2[[1]], "ggplot")
  expect_equal(length(p2), 5) # num_comps + 1
  expect_warning(
    plot_pred(fit, verbose = FALSE),
    "Setting pred=NULL is deprecated."
  )
  expect_warning(
    plot_components(fit, verbose = FALSE, draw = FALSE),
    "Setting pred=NULL is deprecated."
  )
})

test_that("plotting works with reduce = NULL", {
  a <- pred(fit = fit, reduce = NULL, verbose = FALSE)
  p1 <- plot_pred(fit, a, alpha = 0.3) # only mean lines
  expect_s3_class(p1, "ggplot")

  # Color by factor
  cb <- c(NA, "SEX", "SEX", "LOC", "SEX")
  p2 <- plot_components(fit, pred = a, color_by = cb, draw = FALSE)
  p_idx <- sample.int(5, size = 1)
  expect_s3_class(p2[[p_idx]], "ggplot")

  # Expect good errors with invalid input
  expect_error(
    plot_components(fit, pred = a, color_by = c(NA, NA), draw = FALSE),
    "color_by has length 2, but its length should be 5 or one"
  )
  expect_error(
    plot_components(fit, pred = a, color_by = "banana", draw = FALSE),
    "Variable with name 'banana' not found"
  )
})

test_that("plotting works with defaults and given x", {
  x_pred <- new_x(fit, x_values = seq(0, 100, by = 2.5))
  a <- pred(fit = fit, x = x_pred, verbose = FALSE)
  p1 <- plot_pred(fit, pred = a) # nice and smooth plot
  p2 <- plot_components(fit, pred = a, draw = FALSE)
  p_idx <- sample.int(5, size = 1)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2[[p_idx]], "ggplot")
})

test_that("plotting works with defaults and given x, and reduce = NULL", {
  x_pred <- new_x(fit, x_values = seq(0, 100, by = 2.5))
  a <- pred(
    fit = fit, x = x_pred, reduce = NULL, draws = c(2, 5),
    verbose = FALSE
  )
  p1 <- plot_pred(fit, pred = a)
  p2 <- plot_components(fit, pred = a, color_by = "LOC", draw = FALSE)
  p_idx <- sample.int(5, size = 1)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2[[p_idx]], "ggplot")
})