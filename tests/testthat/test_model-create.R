library(lgpr)

# -------------------------------------------------------------------------

context("Creating a model")

test_that("a model can be created using verbose mode", {
  f <- y ~ gp(age) + categ(id) + categ(id) * gp(age)
  expect_output(create_model(f, testdata_001, verbose = TRUE))
})

test_that("created model is valid", {
  f <- y ~ gp(age) + categ(id) + categ(id) * gp(age)
  m <- create_model(f, testdata_001)
  a <- m@model_formula
  expect_equal(.class2(a), "lgpformula")
  expect_equal(length(a@terms@summands), 3)
  expect_equal(length(rhs_variables(a@terms)), 4)
  expect_true(validObject(a))
})

test_that("cannot create nb model where sample_f is FALSE", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  reason <- "sample_f must be TRUE when likelihood is nb"
  expect_error(create_model(
    formula = y ~ gp(age) + categ(sex) * gp(age),
    likelihood = "nb",
    data = dat,
    sample_f = FALSE
  ), reason)
})

test_that("prior can be parsed from stan_input", {
  f <- y ~ gp(age) + categ(id) + categ(id) * gp(age)
  m <- create_model(f, testdata_001)
  df <- prior_to_df(m@stan_input)
  expect_equal(class(df), "data.frame")
})

test_that("quotes can be used in gp(), gp_ns() etc", {
  f <- y ~ gp("age") + categ("id") + categ(id) * gp_ns("age")
  m <- create_model(f, testdata_001)
  a <- m@model_formula
  c <- .class2(a)
  expect_equal(c, "lgpformula")
})

test_that("a heterogeneous component can be added", {
  dat <- testdata_001
  m <- create_model(y ~ het(id) * zs(sex) * gp(age) + categ(id), dat)
  si <- m@stan_input
  idx_heter <- si$components[1, 4]
  expect_equal(idx_heter, 1)
})

test_that("a heterogeneity component must take a categorical covariate", {
  dat <- testdata_001
  expect_error(
    create_model(y ~ het(age) * zs(sex) + categ(id), dat),
    "argument for <het> must be a name of a categorical covariate"
  )
})

test_that("an uncertainty component must take a categorical covariate", {
  dat <- testdata_001
  expect_error(
    create_model(y ~ unc(age) * zs(sex) + categ(id), dat),
    "argument for <unc> must be a name of a categorical covariate"
  )
})

test_that("categ() and zs() must take a categorical covariate", {
  dat <- testdata_001
  expect_error(
    create_model(y ~ gp(age) * zs(sex) + categ(age), dat),
    "argument for <categ> must be a name of a categorical covariate"
  )
  expect_error(
    create_model(y ~ gp(age) * zs(sex) + zs(age), dat),
    "argument for <zs> must be a name of a categorical covariate"
  )
})

test_that("gp(), gp_ns() and gp_vm() must take a continuous covariate", {
  expect_error(
    create_model(y ~ gp(id), testdata_001),
    "argument for <gp> must be a name of a continuous covariate"
  )
  expect_error(
    create_model(y ~ gp_ns(id), testdata_001),
    "argument for <gp_ns> must be a name of a continuous covariate"
  )
  expect_error(
    create_model(y ~ gp_vm(sex), testdata_001),
    "argument for <gp_vm> must be a name of a continuous covariate"
  )
})

test_that("one term can contain at most one expression of each type", {
  expect_error(
    create_model(y ~ gp(id) * gp_ns(age), testdata_001),
    "cannot have more than one gp"
  )
  expect_error(
    create_model(y ~ categ(id) * zs(sex), testdata_001),
    "each term can contain at most one"
  )
  expect_error(
    create_model(y ~ het(id) * gp_ns(age) * het(sex), testdata_001),
    "cannot have more than one"
  )
})

test_that("terms with unc() or het() must have other expressions", {
  expect_error(
    create_model(y ~ het(id), testdata_001),
    "there must be one"
  )
  expect_error(
    create_model(y ~ unc(id), testdata_001),
    "there must be one"
  )
})

test_that("covariate must be same in unc() and het() expression", {
  dat <- testdata_001
  expect_error(
    create_model(y ~ het(id) * unc(sex) * gp(age) + gp(dis_age), dat),
    "Names of the covariates in"
  )
  expect_error(
    create_model(y ~ het(id) * gp(age) + unc(sex) * gp(dis_age), dat),
    "expressions must have the same categorical covariate in every term"
  )
})

test_that("a formula term cannot have more than 4 expressions", {
  dat <- testdata_001
  expect_error(
    create_model(y ~ gp(age) * gp(age) * gp(sex) * het(id) * unc(id), dat),
    "term can have at most four expressions! found = 5"
  )
})

test_that("an lgpmodel has correct list fields for stan input", {
  m <- create_model(
    formula = y ~ gp(age) + zs(id),
    data = testdata_001,
    options = list(delta = 1e-5)
  )
  L <- length(m@stan_input)
  expect_equal(L, 45)
})

test_that("creating an lgpmodel errors with invalid data", {
  expect_error(
    create_model(notvar ~ gp(age) + categ(id), testdata_001),
    "variable 'notvar' not found in <data>"
  )
  expect_error(
    create_model(y ~ gp(notvar) + categ(id), testdata_001),
    "variable 'notvar' not found in <data>"
  )
  expect_error(
    create_model(y ~ gp(age) + categ(id), "notdata"),
    "<data> must be a data.frame"
  )
})

test_that("the num_trials argument works correctly", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  m <- create_model(y ~ gp(age) + zs(id), dat,
    likelihood = "binomial",
    num_trials = 10
  )
  nt <- as.numeric(m@stan_input$y_num_trials)
  expect_equal(nt, rep(10, times = 24))
  reason <- "Invalid length of <num_trials>"
  expect_error(
    create_model(y ~ gp(age) + zs(id), dat,
      likelihood = "binomial",
      num_trials = c(10, 4, 6)
    ),
    reason
  )
})

test_that("only the covariates required by the model go to stan data", {
  m <- create_model(y ~ categ(sex) + zs(id), testdata_001)
  to_stan <- m@stan_input
  expect_equal(to_stan$num_cov_cat, 2)
  expect_equal(dim(to_stan$x_cont), c(0, 24))
})

test_that("covariate types are correctly parsed", {
  m <- create_model(y ~ gp(age) + categ(id) * gp(age) + zs(sex) +
    gp_ns(dis_age), testdata_001)
  to_stan <- m@stan_input
  expect_equal(to_stan$num_cov_cat, 2)
  expect_equal(to_stan$num_cov_cont, 2)
  ts <- as.numeric(to_stan$x_cat_num_levels)
  expect_equal(ts, c(4, 2))
  expect_equal(sum(to_stan$x_cont_mask), 12)
})


test_that("cannot have a continuous covariate with zero variance", {
  newdat <- testdata_001
  newdat$age <- -1
  reason <- "have zero variance"
  expect_error(create_model(y ~ gp(age) + zs(id), newdat), reason)
})

test_that("cannot have NaNs on differing rows", {
  newdat <- testdata_001
  newdat$new_x <- rev(newdat$dis_age)
  reason <- paste0(
    "NaNs of the continuous covariate must be on the same rows.",
    " Found discrepancy between dis_age and new_x"
  )
  expect_error(
    create_model(
      formula = y ~ het(id) * gp(dis_age) + het(id) * gp(new_x),
      data = newdat
    ),
    reason
  )
})

test_that("cannot have missing values for a factor", {
  newdat <- testdata_001
  newdat$sex[3:5] <- NA
  expect_error(
    create_model(formula = y ~ age + sex, data = newdat),
    "missing values for factor 'sex'"
  )
})

test_that("cannot have wrong length of prior list", {
  prior_alpha <- list(log_normal(1, 1), normal(0, 1))
  reason <- "should have length 1 or 3! found = 2"
  expect_error(
    create_model(y ~ id + age + sex,
      data = testdata_001,
      prior = list(alpha = prior_alpha)
    ),
    reason
  )
})

test_that("verbose mode prints output", {
  expect_output(
    create_model(y ~ id + age + sex,
      data = testdata_001,
      verbose = TRUE
    )
  )
})
