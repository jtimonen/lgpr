library(lgpr)

# Create test data
age <- c(10, 20, 30, 10, 20, 30)
id <- as.factor(c(1, 1, 1, 2, 2, 2))
sex <- as.factor(c("Male", "Male", "Male", "Female", "Female", "Female"))
dis_age <- c(12 - age[1:3], NA, NA, NA)
y <- c(1, 2, 4, 10, 0, 4)
dat <- data.frame(id, age, dis_age, sex, y)

# -------------------------------------------------------------------------

context("Input parsing functions")

test_that("parse_formula works with a single lgpexpr", {
  f <- y ~ gp(x)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("response variable cannot be a covariate", {
  f <- y ~ gp(x) + categ(y) * gp(t)
  msg <- "the response variable cannot be also a covariate"
  expect_error(parse_formula(f), msg)
})

test_that("parse_formula works with a single lgpterm", {
  f <- y ~ gp(x) * categ(z)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("two lgpterms can be summed", {
  f <- y ~ gp(x) * categ(z) + gp(x) * zerosum(z)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpterm and lgpexpr can be summed", {
  f <- y ~ gp(x) * zerosum(z) + categ("z")
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpexpr and lgpterm can be summed", {
  f <- y ~ gp(x) + categ("z") * gp_warp(aa)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("parse_formula throws error when input is not a formula", {
  reason <- "must have class formula"
  expect_error(parse_formula("a + b "), reason)
})

test_that("parse_formula throws error if invalid function or covariate", {
  expect_error(parse_formula(y ~ notafunction(x)), "<fun> must be one of")
  expect_error(parse_formula(y ~ gp("")), "covariate name cannot be empty")
  expect_error(
    parse_formula(y ~ gp(x) * gp(y) * gp(z)),
    "the response variable cannot be also a covariate"
  )
  expect_error(
    parse_formula(y ~ x + a),
    "expression must contain exactly one opening and closing parenthesis"
  )
  expect_error(
    parse_formula(y ~ gp(x)(y)),
    "expression must contain exactly one opening and closing parenthesis"
  )
})

test_that("parse_options does not need arguments", {
  a <- parse_options()
  e <- c("is_generated_skipped", "is_f_sampled", "delta")
  expect_equal(names(a), e)
})

test_that("parse_likelihood can be used", {
  list_y <- list(y_cont = c(1, 2))
  a <- parse_likelihood("gaussian", NULL, NULL, list_y)
  e <- c("obs_model", "y_num_trials", "c_hat", "is_likelihood_skipped")
  expect_equal(names(a), e)
})

test_that("parse_likelihood errors correctly", {
  list_y <- list(y_cont = c(1, 2))
  expect_error(parse_likelihood("gaussian", 0, NULL, list_y))
  expect_error(parse_likelihood("gaussian", NULL, 0, list_y))
  list_y <- list(y_disc = c(1, 2))
  expect_error(parse_likelihood("nb", c(1, 2, 3), NULL, list_y))
})

test_that("parse_likelihood works correctly with binomial likelihood", {
  list_y <- list(y_disc = c(1, 2))
  a <- parse_likelihood("binomial", 1.2, 9, list_y)
  b <- parse_likelihood("binomial", NULL, 9, list_y)
  expect_equal(a$obs_model, 4)
  expect_equal(length(a$y_num_trials), 2)
  diff <- abs(b$c_hat - (-1.609438))
  expect_lt(max(diff), 1e-6)
  expect_error(parse_likelihood("binomial", NULL, c(1, 1, 1), c(1, 2)))
  expect_error(parse_likelihood("binomial", c(1, 1, 1), NULL, c(1, 2)))
})

test_that("y_scaling is created and and applied", {
  f <- parse_formula(y ~ gp(age) + zerosum(id))
  parsed <- parse_response(dat, "gaussian", f)
  yts <- parsed$to_stan
  expect_equal(names(parsed), c("to_stan", "scaling"))
  expect_equal(yts$num_obs, 6)
  expect_true(class(parsed$scaling) == "lgpscaling")

  # Check zero mean and unit variance of y_cont
  y <- as.numeric(yts$y_cont)
  d1 <- abs(mean(y) - 0)
  d2 <- abs(stats::sd(y) - 1)
  expect_lt(d1, 1e-6)
  expect_lt(d2, 1e-6)
})

test_that("cannot have negative response with NB observation model", {
  newdat <- dat
  newdat$y <- c(-1, -9, 3, 2, 4, 1)
  f <- parse_formula(y ~ gp(age) + zerosum(id))
  reason <- "cannot be negative with this observation model"
  expect_error(parse_response(newdat, "nb", f), reason)
})

test_that("cannot have a response with zero variance (gaussian obs model)", {
  newdat <- dat
  newdat$y <- c(1, 1, 1, 1, 1, 1)
  f <- parse_formula(y ~ gp(age) + zerosum(id))
  reason <- "have zero variance"
  expect_error(parse_response(newdat, "gaussian", f), reason)
})
