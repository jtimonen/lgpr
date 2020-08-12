library(lgpr)


# -------------------------------------------------------------------------

context("Input parsing: parsing formulas")

test_that("parse_formula returns correct object", {
  f <- y ~ gp(x) + categ(u) + categ(z) * gp(x) + gp_ns(age) + zerosum(q)
  a <- parse_formula(f)
  expect_equal(.class2(a), "lgpformula")
  expect_equal(length(a@terms@summands), 5)
  expect_equal(length(rhs_variables(a@terms)), 6)
  expect_true(validObject(a))
  expect_gt(nchar(as.character(a)), 10)
})

test_that("quotes can be used in gp(), gp_ns() etc", {
  f <- y ~ gp("x") + categ("u") + categ("z") * gp("x") +
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
  f <- y ~ gp(x) + categ("z") * gp_ns(aa)
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

test_that("parse_response creates y_scaling and applies it", {
  f <- parse_formula(y ~ gp(age) + zerosum(id))
  parsed <- parse_response(dat, "gaussian", f)
  yts <- parsed$to_stan
  expect_equal(names(parsed), c("to_stan", "scaling"))
  expect_equal(yts$num_obs, 6)
  expect_true(class(parsed$scaling) == "lgpscaling")

  # Check zero mean and unit variance of y_cont
  y <- as.numeric(yts$y_cont_norm)
  d1 <- abs(mean(y) - 0)
  d2 <- abs(stats::sd(y) - 1)
  expect_lt(d1, 1e-6)
  expect_lt(d2, 1e-6)
})

# -------------------------------------------------------------------------

context("Input parsing: creating an lgpmodel")

test_that("an lgpmodel can be created", {
  m <- lgp_model(y ~ gp(age) + zerosum(id), dat)
  expect_true(class(m) == "lgpmodel")
})

# test_that("an lgpmodel has correct number of list fields for stan input", {
#  m <- lgp_model(y ~ gp(age) + mask(id), dat, options = list(delta = 1e-5))
#  L <- length(names(m@stan_input))
#  expect_equal(L, 22)
#  expect_equal(m@stan_input$delta, 1e-5)
# })

test_that("createing lgpmodel errors correctly", {
  expect_error(lgp_model(notvar ~ gp(age) + categ(id), dat))
  expect_error(lgp_model(y ~ gp(notvar) + categ(id), dat))
  expect_error(lgp_model(y ~ gp(age) + categ(id), "notdata"))
  expect_error(lgp_model(y ~ gp(age) + categ(id), dat, likelihood = "binomial"))
})

# Create larger test data
age <- c(10, 20, 30, 10, 20, 30)
id <- as.factor(c(1, 1, 1, 2, 2, 2))
sex <- as.factor(c("Male", "Male", "Male", "Female", "Female", "Female"))
dis_age <- c(12 - age[1:3], NA, NA, NA)
y <- c(1, 2, 4, 10, 0, 4)
dat <- data.frame(id, age, dis_age, sex, y)

test_that("only the covariates required by the model go to stan data", {
  m <- lgp_model(y ~ gp(sex) + zerosum(id), dat)
  to_stan <- m@stan_input
  expect_equal(to_stan$num_cov_cat, 2)
  expect_equal(dim(to_stan$x_cont), c(0, 6))
})

test_that("covariate types are correctly parsed", {
  m <- lgp_model(y ~ gp(age) + categ(id) * gp(age) + zerosum(sex) +
    gp_ns(dis_age), dat)
  to_stan <- m@stan_input
  expect_equal(to_stan$num_cov_cat, 2)
  expect_equal(to_stan$x_cat_num_levels, c(2, 2))
  expect_equal(to_stan$num_cov_cont, 2)
  expect_equal(dim(to_stan$x_cont_normalized), c(2, 6))
  expect_equal(sum(to_stan$x_cont_mask), 3)
})
