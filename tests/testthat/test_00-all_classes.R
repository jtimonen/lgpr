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
