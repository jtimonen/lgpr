library(lgpr)
source("helpers/SW.R")

# -------------------------------------------------------------------------

context("Validation of S4 class objects")

test_that("lgpexpr validation works correctly", {
  a <- lgpexpr(fun = "gp", covariate = "x")
  expect_true(validate_lgpexpr(a))
  b <- a
  b@fun <- "moi"
  msg <- validate_lgpexpr(b)
  expect_error(stop(msg), "<fun> must be one of")
  c <- a
  c@covariate <- ""
  msg <- validate_lgpexpr(c)
  expect_error(stop(msg), "covariate name cannot be empty")
})

test_that("lgpformula validation works correctly", {
  a <- create_model.formula(as.formula("y ~ gp(x) + zs(a)"), NULL)
  expect_true(validate_lgpformula(a))
  b <- a
  b@y_name <- "x"
  msg <- validate_lgpformula(b)
  expect_error(stop(msg), "response variable cannot be also")
})

test_that("lgpscaling and its validation work correctly", {
  # Test validation
  x <- stats::rnorm(1000, mean = -300, sd = 2)
  x[1:500] <- NA
  a <- create_scaling(x, "test_var")
  expect_true(validate_lgpscaling(a))
  b <- a
  b@var_name <- ""
  msg <- validate_lgpscaling(b)
  expect_error(stop(msg), "variable name length must be at least 1")

  # Test that scaling works correctly
  y <- apply_scaling(a, x)
  my <- mean(y, na.rm = TRUE)
  sy <- stats::sd(y, na.rm = TRUE)
  expect_equal(my, 0.0)
  expect_equal(sy, 1.0)

  # Test that inverse scaling works correctly
  x_rec <- apply_scaling(a, y, inverse = TRUE)
  max_err <- max(abs(x - x_rec), na.rm = TRUE)
  expect_lt(max_err, 1e-6)
})

test_that("lgpfit and GaussianPrediction validation works correctly", {
  SW({
    fit <- example_fit(refresh = 0, quiet = TRUE)
  })
  expect_true(validate_lgpfit(fit))
  a <- pred(fit, verbose = FALSE)
  expect_true(validate_GaussianPrediction(a))
})

test_that("lgpfit Prediction validation works correctly", {
  SW({
    fit <- example_fit(refresh = 0, quiet = TRUE, likelihood = "nb")
  })
  expect_true(validate_lgpfit(fit))
  a <- pred(fit, verbose = FALSE, reduce = NULL)
  expect_true(validate_Prediction(a))
})
