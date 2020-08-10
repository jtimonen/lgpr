library(lgpr)


# -------------------------------------------------------------------------

context("Input parsing: parsing formulas")

test_that("parse_formula returns correct object", {
  f <- y ~ gp(x) + categorical(u) + mask(z) * gp(x) + gp_ns(age) + zerosum(q)
  a <- parse_formula(f)
  expect_equal(.class2(a), "lgpformula")
  expect_equal(length(a@terms@summands), 5)
  expect_equal(length(rhs_variables(a@terms)), 6)
  expect_true(validObject(a))
  expect_gt(nchar(as.character(a)), 10)
})

test_that("quotes can be used in gp(), gp_ns() etc", {
  f <- y ~ gp("x") + categorical("u") + mask("z") * gp("x") +
    gp_ns("age") + zerosum("q")
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("parse_formula works with a single lgpexpr", {
  f <- y ~ gp(x)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("response variable cannot be a covariate", {
  f <- y ~ gp(x) + mask(y) * gp(t)
  msg <- "the response variable cannot be also a covariate"
  expect_error(parse_formula(f), msg)
})

test_that("parse_formula works with a single lgpterm", {
  f <- y ~ gp(x) * mask(z)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("two lgpterms can be summed", {
  f <- y ~ gp(x) * mask(z) + gp(x) * mask(z)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpterm and lgpexpr can be summed", {
  f <- y ~ gp(x) * mask(z) + mask("z")
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpexpr and lgpterm can be summed", {
  f <- y ~ gp(x) + mask("z") * gp_ns(aa)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("parse_formula throws error when input is not a formula", {
  expect_error(parse_formula("a + b "))
})

test_that("parse_formula throws error if invalid function or covariate", {
  expect_error(parse_formula(y ~ notafunction(x)))
  expect_error(parse_formula(y ~ gp("")))
  expect_error(parse_formula(y ~ gp(x) * gp(y) * gp(z)))
  expect_error(parse_formula(y ~ x + a))
  expect_error(parse_formula(y ~ gp(x)(y)))
})

# -------------------------------------------------------------------------

context("Input parsing: parsing options")

test_that("parse_options does not need arguments", {
  a <- parse_options()
  e <- c("is_verbose", "is_generated_skipped", "is_f_sampled", "delta")
  expect_equal(names(a), e)
})

test_that("parse_disease_options does not need arguments", {
  a <- parse_disease_options()
  e <- c("is_uncrt", "is_heter", "is_vm_used", "vm_params")
  expect_equal(names(a), e)
})

test_that("vm_params gets correct dimensions", {
  a <- parse_disease_options()
  b <- parse_disease_options(list(uncertain = TRUE, vm_params = NA))
  expect_equal(dim(a$vm_params), c(1, 2))
  expect_equal(dim(b$vm_params), c(0, 2))
})

# -------------------------------------------------------------------------

context("Input parsing: parsing likelihood")

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

# -------------------------------------------------------------------------

context("Input parsing: parsing data")

# Create test data
age <- c(10, 20, 30, 10, 20, 30)
id <- c(1, 1, 1, 2, 2, 2)
y <- c(9.3, 1.2, 3.2, 2.3, 4.1, 1)
dat <- data.frame(age, id, y)

test_that("parse_response creates y_scaling", {
  f <- parse_formula(y ~ gp(age) + mask(id))
  parsed <- parse_response(dat, "gaussian", f)
  expect_equal(names(parsed), c("y_to_stan", "y_scaling"))
  expect_equal(parsed$y_to_stan$num_obs, 6)
  expect_true(class(parsed$y_scaling) == "lgpscaling")
})

# -------------------------------------------------------------------------

context("Input parsing: creating an lgpmodel")

test_that("an lgpmodel can be created", {
  m <- lgp_model(y ~ gp(age) + mask(id), dat)
  expect_true(class(m) == "lgpmodel")
})

test_that("an lgpmodel has correct number of list fields for stan input", {
  m <- lgp_model(y ~ gp(age) + mask(id), dat, options = list(delta = 1e-5))
  L <- length(names(m@stan_input))
  expect_equal(L, 15)
  expect_equal(m@stan_input$delta, 1e-5)
})

test_that("createing lgpmodel errors correctly", {
  expect_error(lgp_model(notvar ~ gp(age) + mask(id), dat))
  expect_error(lgp_model(y ~ gp(age) + mask(id), "notdata"))
  expect_error(lgp_model(y ~ gp(age) + mask(id), dat, likelihood = "binomial"))
})
