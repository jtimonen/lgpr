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

test_that("parse_formula translated simple formulas", {
  f <- parse_formula(y ~ age + sex, dat)
  expect_equal(as.character(f), "y ~ gp(age) + zs(sex)")
  f <- parse_formula(y ~ age + age | sex, dat)
  expect_equal(as.character(f), "y ~ gp(age) + gp(age) * zs(sex)")
  expect_output(show(f@terms))
  expect_output(show(f@terms@summands[[1]]))
  reason <- "variable 'notvar' not found in <data>"
  expect_error(parse_formula(y ~ age + notvar, dat), reason)
})

test_that("parse_formula_advanced works with a single lgpexpr", {
  f <- y ~ gp(x)
  c <- .class2(parse_formula_advanced(f))
  expect_equal(c, "lgpformula")
})

test_that("response variable cannot be a covariate", {
  f <- y ~ gp(x) + categ(y) * gp(t)
  msg <- "the response variable cannot be also a covariate"
  expect_error(parse_formula_advanced(f), msg)
})

test_that("parse_formula_advanced works with a single lgpterm", {
  f <- y ~ gp(x) * categ(z)
  c <- .class2(parse_formula_advanced(f))
  expect_equal(c, "lgpformula")
})

test_that("two lgpterms can be summed", {
  f <- y ~ gp(x) * categ(z) + gp(x) * zs(z)
  c <- .class2(parse_formula_advanced(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpterm and lgpexpr can be summed", {
  f <- y ~ gp(x) * zs(z) + categ("z")
  c <- .class2(parse_formula_advanced(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpexpr and lgpterm can be summed", {
  f <- y ~ gp(x) + categ("z") * gp_ns(aa)
  c <- .class2(parse_formula_advanced(f))
  expect_equal(c, "lgpformula")
})

test_that("parse_formula_advanced throws error when input is not a formula", {
  reason <- "Allowed classes are"
  expect_error(parse_formula_advanced("a + b "), reason)
})

test_that("parse_formula_advanced errors with invalid function or covariate", {
  expect_error(parse_formula_advanced(y ~ notfun(x)), "<fun> must be one of")
  expect_error(parse_formula_advanced(y ~ gp("")), "cannot be empty")
  expect_error(
    parse_formula_advanced(y ~ gp(x) * gp(y) * gp(z)),
    "the response variable cannot be also a covariate"
  )
  expect_error(
    parse_formula_advanced(y ~ gp(x) + gp((a))),
    "expression must contain exactly one opening and closing parenthesis"
  )
  expect_error(
    parse_formula_advanced(y ~ gp(x)(y)),
    "expression must contain exactly one opening and closing parenthesis"
  )
})


test_that("parse_formula_advanced throws error if mixing syntaxes", {
  reason <- "Be sure not to mix the advanced and simple formula syntaxes"
  expect_error(parse_formula_advanced(5 ~ koira + gp(r)), reason)
})


test_that("parse_options does not need arguments", {
  a <- parse_options()
  e <- c("is_generated_skipped", "delta")
  expect_equal(names(a), e)
})

test_that("parse_likelihood can be used", {
  list_y <- list(y_cont = c(1, 2))
  a <- parse_likelihood("gaussian", NULL, NULL, list_y, TRUE)
  e <- c("obs_model", "y_num_trials", "c_hat", "is_f_sampled")
  expect_equal(names(a), e)
})

test_that("parse_likelihood errors correctly", {
  list_y <- list(y_cont = c(1, 2))
  reason <- "Only give the c_hat argument if observation model is not Gaussian"
  expect_error(parse_likelihood("gaussian", 0, NULL, list_y, FALSE), reason)
  reason <- "<num_trials> argument if likelihood is binomial or beta-binomial"
  expect_error(parse_likelihood("gaussian", NULL, 0, list_y, FALSE), reason)
  list_y <- list(y_disc = c(1, 2))
  reason <- "Invalid length of <c_hat>"
  expect_error(parse_likelihood("nb", c(1, 2, 3), NULL, list_y, FALSE), reason)
})

test_that("parse_likelihood works correctly with binomial likelihood", {
  list_y <- list(y_disc = c(1, 2))
  a <- parse_likelihood("binomial", 1.2, 9, list_y, TRUE)
  b <- parse_likelihood("binomial", NULL, 9, list_y, TRUE)
  expect_equal(a$obs_model, 4)
  expect_equal(length(a$y_num_trials), 2)
  diff <- abs(b$c_hat - (-1.609438))
  expect_lt(max(diff), 1e-6)
  reason <- "Invalid length of <num_trials>"
  expect_error(
    parse_likelihood("binomial", NULL, c(1, 1, 1), list_y, FALSE),
    reason
  )
  reason <- "Invalid length of <c_hat>"
  expect_error(
    parse_likelihood("binomial", c(1, 1, 1), NULL, list_y, FALSE),
    reason
  )
})

test_that("y_scaling is created and and applied", {
  f <- parse_formula_advanced(y ~ gp(age) + zs(id))
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
  f <- parse_formula_advanced(y ~ gp(age) + zs(id))
  reason <- "<y> must have only non-negative values"
  expect_error(parse_response(newdat, "nb", f), reason)
})

test_that("cannot have a response with zero variance (gaussian obs model)", {
  newdat <- dat
  newdat$y <- c(1, 1, 1, 1, 1, 1)
  f <- parse_formula_advanced(y ~ gp(age) + zs(id))
  reason <- "have zero variance"
  expect_error(parse_response(newdat, "gaussian", f), reason)
})
