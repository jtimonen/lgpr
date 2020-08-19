library(lgpr)

# Create test data
age <- c(10, 20, 30, 10, 20, 30)
id <- as.factor(c(1, 1, 1, 2, 2, 2))
sex <- as.factor(c("Male", "Male", "Male", "Female", "Female", "Female"))
dis_age <- c(12 - age[1:3], NA, NA, NA)
y <- c(1, 2, 4, 10, 0, 4)
dat <- data.frame(id, age, dis_age, sex, y)

# -------------------------------------------------------------------------

context("Input parsing")

test_that("lgp_model can be used", {
  f <- y ~ gp(age) + categ(id) + categ(id) * gp(age)
  m <- lgp_model(f, dat)
  a <- m@model_formula
  expect_equal(.class2(a), "lgpformula")
  expect_equal(length(a@terms@summands), 3)
  expect_equal(length(rhs_variables(a@terms)), 4)
  expect_true(validObject(a))
})

test_that("lgpmodel and lgpformula have character representations", {
  f <- y ~ gp(age) + categ(id) + categ(id) * gp(age)
  m <- lgp_model(f, dat)
  a <- m@model_formula
  expect_gt(nchar(as.character(m)), 10)
  expect_gt(nchar(as.character(a)), 10)
})

test_that("quotes can be used in gp(), gp_warp() etc", {
  f <- y ~ gp("age") + categ("id") + categ(id) * gp_warp("age")
  m <- lgp_model(f, dat)
  a <- m@model_formula
  c <- .class2(a)
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
  e <- c("is_verbose", "is_generated_skipped", "is_f_sampled", "delta")
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

test_that("a heterogeneous component can be added", {
  m <- lgp_model(y ~ heter(id) * zerosum(sex) * gp(age) + categ(id), dat)
  si <- m@stan_input
  idx_heter <- si$components[1, 4]
  expect_equal(idx_heter, 1)
})

test_that("an uncertain component can be added", {
  m <- lgp_model(y ~ uncrt(id) * zerosum(sex) * gp(age) + categ(id), dat)
  si <- m@stan_input
  idx_heter <- si$components[1, 4]
  idx_uncrt <- si$components[1, 7]
  expect_equal(idx_heter, 0)
  expect_equal(idx_uncrt, 1)
})

test_that("a heterogeneity component must take a categorical covariate", {
  expect_error(
    lgp_model(y ~ heter(age) * zerosum(sex) + categ(id), dat),
    "argument for <heter> must be a name of a categorical covariate"
  )
})

test_that("an uncertainty component must take a categorical covariate", {
  expect_error(
    lgp_model(y ~ uncrt(age) * zerosum(sex) + categ(id), dat),
    "argument for <uncrt> must be a name of a categorical covariate"
  )
})

test_that("categ() and zerosum() must take a categorical covariate", {
  expect_error(
    lgp_model(y ~ gp(age) * zerosum(sex) + categ(age), dat),
    "argument for <categ> must be a name of a categorical covariate"
  )
  expect_error(
    lgp_model(y ~ gp(age) * zerosum(sex) + zerosum(age), dat),
    "argument for <zerosum> must be a name of a categorical covariate"
  )
})

test_that("gp(), gp_warp() and gp_warp_vm must take a continuous covariate", {
  expect_error(
    lgp_model(y ~ gp(id), dat),
    "argument for <gp> must be a name of a continuous covariate"
  )
  expect_error(
    lgp_model(y ~ gp_warp(id), dat),
    "argument for <gp_warp> must be a name of a continuous covariate"
  )
  expect_error(
    lgp_model(y ~ gp_warp_vm(sex), dat),
    "argument for <gp_warp_vm> must be a name of a continuous covariate"
  )
})

test_that("one term can contain at most one expression of each type", {
  expect_error(
    lgp_model(y ~ gp(id) * gp_warp(age), dat),
    "cannot have more than one gp"
  )
  expect_error(
    lgp_model(y ~ categ(id) * zerosum(sex), dat),
    "each term can contain at most one"
  )
  expect_error(
    lgp_model(y ~ heter(id) * gp_warp(age) * heter(sex), dat),
    "cannot have more than one 'heter' expression in one term"
  )
})

test_that("terms with uncrt() or heter() must have other expressions", {
  expect_error(
    lgp_model(y ~ heter(id), dat),
    "there must be one <gp>, <gp_warp> or <gp_warp_vm> expression"
  )
  expect_error(
    lgp_model(y ~ uncrt(id), dat),
    "there must be one <gp>, <gp_warp> or <gp_warp_vm> expression"
  )
})

test_that("covariate must be same in uncrt() and heter() expression", {
  expect_error(
    lgp_model(y ~ heter(id) * uncrt(sex) * gp(age) + gp(dis_age), dat),
    "expressions must match! Found = "
  )
})

test_that("a formula  term cannot have more than 4 expressions", {
  expect_error(
    lgp_model(y ~ gp(age) * gp(age) * gp(sex) * heter(id) * uncrt(id), dat),
    "term can have at most four expressions! found = 5"
  )
})

# test_that("terms with uncrt() or heter() expression must contain a gp expr", {
#  print("TODO")
# })

test_that("an lgpmodel has correct list fields for stan input", {
  m <- lgp_model(y ~ gp(age) + zerosum(id), dat, options = list(delta = 1e-5))
  found_fields <- sort(names(m@stan_input))
  expected_fields <- sort(stan_list_names())
  for (field in expected_fields) {
    expect_true(!!field %in% found_fields)
  }
  expect_equal(m@stan_input$delta, 1e-5)
})

test_that("creating an lgpmodel errors with invalid data", {
  expect_error(
    lgp_model(notvar ~ gp(age) + categ(id), dat),
    "variable 'notvar' not found in <data>"
  )
  expect_error(
    lgp_model(y ~ gp(notvar) + categ(id), dat),
    "variable 'notvar' not found in <data>"
  )
  expect_error(
    lgp_model(y ~ gp(age) + categ(id), "notdata"),
    "<data> must be a data.frame"
  )
})

test_that("the num_trials argument works correctly", {
  dat$y <- c(1, 3, 2, 7, 3, 1)
  m <- lgp_model(y ~ gp(age) + zerosum(id), dat,
    likelihood = "binomial",
    num_trials = 10
  )
  nt <- as.numeric(m@stan_input$y_num_trials)
  expect_equal(nt, rep(10, times = 6))
  reason <- "Invalid length of <num_trials>"
  expect_error(
    lgp_model(y ~ gp(age) + zerosum(id), dat,
      likelihood = "binomial",
      num_trials = c(10, 4, 6)
    ),
    reason
  )
})

test_that("only the covariates required by the model go to stan data", {
  m <- lgp_model(y ~ categ(sex) + zerosum(id), dat)
  to_stan <- m@stan_input
  expect_equal(to_stan$num_cov_cat, 2)
  expect_equal(dim(to_stan$x_cont), c(0, 6))
})

test_that("covariate types are correctly parsed", {
  m <- lgp_model(y ~ gp(age) + categ(id) * gp(age) + zerosum(sex) +
    gp_warp(dis_age), dat)
  to_stan <- m@stan_input
  expect_equal(to_stan$num_cov_cat, 2)
  expect_equal(to_stan$num_cov_cont, 2)
  ts <- as.numeric(to_stan$x_cat_num_levels)
  expect_equal(ts, c(2, 2))
  expect_equal(sum(to_stan$x_cont_mask), 3)
})


test_that("cannot have a continuous covariate with zero variance", {
  newdat <- dat
  newdat$age <- -1
  reason <- "have zero variance"
  expect_error(lgp_model(y ~ gp(age) + zerosum(id), newdat), reason)
})
