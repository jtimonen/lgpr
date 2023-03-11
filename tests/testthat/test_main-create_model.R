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
  expect_output(show(a))
  expect_equal(.class2(a), "lgpformula")
  expect_equal(length(a@terms@summands), 3)
  expect_equal(length(rhs_variables(a@terms)), 4)
  expect_true(validObject(a))
})

test_that("cannot create nb model where sample_f is FALSE", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  reason <- "<sample_f> must be TRUE when <likelihood> is nb"
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
  m <- create_model(f, testdata_001, prior = list(wrp = igam(14, 5)))
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
    "Cannot have more than one gp"
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

test_that("if using z|x in formula it is informed that they are wrong way", {
  msg <- paste0(
    "If you used the | syntax in your model formula, make sure ",
    "that the covariate on the left side of | is continuous ",
    "and the one on the right side is categorical."
  )
  expect_error(
    create_model(y ~ age + id | age, testdata_001),
    msg
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
    "must be an object of type <data.frame>"
  )
})

test_that("the num_trials argument works correctly", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  m <- create_model(y ~ gp(age) + zs(id), dat,
    likelihood = "binomial",
    num_trials = 10
  )
  si <- get_stan_input(m)
  nt <- as.numeric(si$y_num_trials)
  expect_equal(nt, rep(10, times = 24))
})

test_that("only the covariates required by the model go to stan data", {
  m <- create_model(y ~ categ(sex) + zs(id), testdata_001)
  si <- get_stan_input(m)
  expect_equal(si$num_cov_cat, 2)
  expect_equal(dim(si$x_cont), c(0, 24))
})

test_that("covariate types are correctly parsed", {
  m <- create_model(y ~ gp(age) + categ(id) * gp(age) + zs(sex) +
    gp_ns(dis_age), testdata_001, prior = list(wrp = igam(14, 5)))
  si <- get_stan_input(m)
  expect_equal(si$num_cov_cat, 2)
  expect_equal(si$num_cov_cont, 2)
  ts <- as.numeric(si$x_cat_num_levels)
  expect_equal(ts, c(4, 2))
  expect_equal(sum(si$x_cont_mask), 12)
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

test_that("warning is given about default prior for wrp", {
  m <- create_model(y ~ gp_vm(age) + gp_vm(blood),
    prior = list(wrp = igam(14, 5)),
    testdata_001
  )
  expect_warning(
    create_model(y ~ gp_vm(age), testdata_001)
  )
})


test_that("simple formulas are translated to advanced", {
  m <- create_model(y ~ age + sex, testdata_001)
  expect_equal(as.character(m@model_formula), "y ~ gp(age) + zs(sex)")
  m <- create_model(y ~ age + age | sex, testdata_001)
  form <- m@model_formula
  expect_equal(as.character(form), "y ~ gp(age) + gp(age) * zs(sex)")
  expect_output(show(form@terms))
  expect_output(show(form@terms@summands[[1]]))
  reason <- "variable 'notvar' not found in <data>"
  expect_error(create_model(y ~ age + notvar, testdata_001), reason)
})

test_that("response variable cannot be also a covariate", {
  f <- y ~ gp(x) + categ(y) * gp(t)
  msg <- "the response variable cannot be also a covariate"
  expect_error(create_model(f, testdata_001), msg)
})

test_that("advanced formula parsing works with a single lgpexpr", {
  m <- create_model(y ~ gp(age), testdata_001)
  types <- as.numeric(component_info(m)$type)
  cn <- component_names(m)
  expect_equal(cn, c("gp(age)"))
  expect_equal(types, 1)
})

test_that("two lgpterms can be summed", {
  m <- create_model(age ~ gp(y) * categ(sex) + gp(y) * zs(id), testdata_001)
  types <- as.numeric(component_info(m)$type)
  cn <- component_names(m)
  expect_equal(cn, c("gp(y)*categ(sex)", "gp(y)*zs(id)"))
  expect_equal(types, c(2, 2))
})

test_that("lgpterm and lgpexpr can be summed", {
  m <- create_model(y ~ gp(age) * zs(sex) + categ("id"), testdata_001)
  types <- as.numeric(component_info(m)$type)
  cn <- component_names(m)
  expect_equal(cn, c("gp(age)*zs(sex)", "categ(id)"))
  expect_equal(types, c(2, 0))
})

test_that("lgpexpr and lgpterm can be summed", {
  m <- create_model(y ~ gp(age) + categ(sex) * gp_ns(dis_age), testdata_001,
    prior = list(wrp = normal(1, 0.1))
  )
  types <- as.numeric(component_info(m)$type)
  ns <- as.numeric(component_info(m)$ns)
  cn <- component_names(m)
  expect_equal(cn, c("gp(age)", "categ(sex)*gp_ns(dis_age)"))
  expect_equal(types, c(1, 2))
  expect_equal(ns, c(0, 1))
})

test_that("error is thrown when formula is not a formula", {
  reason <- "must be an object of type <formula>! Found = <character>"
  expect_error(create_model("a + b ", testdata_001), reason)
})

test_that("formula parsing errors with invalid formula", {
  dat <- testdata_001
  expect_error(create_model(y ~ notfun(age), dat), "<fun> must be one of")
  expect_error(create_model(y ~ gp(""), dat), "covariate name cannot be empty")
  expect_error(
    create_model(y ~ gp(x) + gp((a)), dat),
    "expression must contain exactly one opening and closing parenthesis"
  )
  expect_error(
    create_model(y ~ gp(x)(y), dat),
    "expression must contain exactly one opening and closing parenthesis"
  )
})

test_that("formula parsing throws error if mixing syntaxes", {
  reason <- "Be sure not to mix the advanced and simple formula syntaxes"
  expect_error(create_model(5 ~ koira + gp(r), testdata_001), reason)
})

test_that("parsing likelihood and response errors correctly", {
  dat <- testdata_001
  reason <- "<c_hat> must be NULL when <sample_f> is FALSE"
  expect_error(create_model(y ~ age, dat, c_hat = 0), reason)
  reason <- "<num_trials> must be NULL when <sample_f> is FALSE"
  expect_error(create_model(y ~ age, dat, num_trials = 100), reason)

  dat$y <- round(exp(dat$y))
  reason <- "<num_trials> argument if likelihood is 'binomial' or 'bb'"
  expect_error(create_model(y ~ age, dat,
    num_trials = 0,
    likelihood = "nb"
  ), reason)
  expect_error(
    create_model(y ~ age, likelihood = "Poisson", dat, c_hat = c(1, 2)),
    "Invalid length of <c_hat>"
  )
  expect_error(
    create_model(y ~ age, likelihood = "binomial", dat, num_trials = c(1, 2)),
    "Invalid length of <num_trials>"
  )
  expect_error(
    create_model(y ~ age, likelihood = "binomial", dat, num_trials = c(3)),
    "value of <response> is larger than value of <num_trials> at"
  )
})


test_that("y_scaling is created and and applied, original data staying as is", {
  dat <- testdata_001
  ADD <- 300
  dat$y <- 100 * dat$y + ADD
  m <- create_model(y ~ gp(age) + zs(id), dat)
  y_scl <- m@var_scalings$y
  expect_true(class(y_scl) == "lgpscaling")
  expect_equal(mean(get_data(m)$y), ADD) # should be original

  # Check zero mean and unit variance of y_norm
  y_norm <- get_stan_input(m)$y_norm
  expect_equal(mean(y_norm), 0.0)
  expect_equal(var(y_norm), 1.0)
})

test_that("cannot have negative or non-integer response with NB etc.", {
  dat <- testdata_001
  dat$y <- exp(dat$y)
  expect_error(
    create_model(y ~ age, dat, likelihood = "NB"),
    "<response> must have only integer values"
  )
  dat <- testdata_001
  dat$y <- round(100 * dat$y)
  expect_error(
    create_model(y ~ age, dat, likelihood = "NB"),
    "<response> must have only non-negative values"
  )
})

test_that("a data variable that will be normalized can't have zero variance", {
  # Response
  dat <- testdata_001
  dat$y <- rep(10.3, nrow(dat))
  reason <- "the variable <y> has zero variance"
  expect_error(create_model(y ~ gp(age) + zs(id), dat), reason)

  # Continuous covariate
  dat <- testdata_001
  dat$age <- rep(-3.2, nrow(dat))
  reason <- "the variable <age> has zero variance"
  expect_error(create_model(y ~ gp(age) + zs(id), dat), reason)
})

test_that("a model with uncertain disease age needs prior specified", {
  formula <- y ~ gp(age) + unc(id) * gp_vm(dis_age)
  data <- testdata_001
  reason <- "you must specify 'effect_time_info' in"
  expect_error(create_model(formula = formula, data = data), reason)
})

test_that("can't have only partly missing vals for group when using het()", {
  formula <- y ~ gp(age) + het(id) * gp(dis_age)
  dat <- testdata_001
  dat$dis_age[20] <- 1.1
  reason <- "inconsistent x_cont_mask values for observations where id = 4"
  expect_error(create_model(formula = formula, data = dat), reason)
})

test_that("options can be specified", {
  form <- y ~ age + id
  dat <- testdata_001
  DVAL <- stats::runif(1, 0, 0.01)
  model <- create_model(form, dat, options = list(delta = DVAL))
  expect_equal(model@stan_input$delta, DVAL)
  model <- create_model(form, dat, options = list(vm_params = c(0.1, 0.3)))
  expect_equal(model@stan_input$vm_params, c(0.1, 0.3))
})

# -------------------------------------------------------------------------

context("Defining priors")


test_that("uniform prior is parsed correctly", {
  suppressWarnings({
    p <- uniform()
    num <- prior_to_num(p)
    expect_equal(num$prior, c(1, 0))
  })
})

test_that("normal prior is parsed correctly", {
  p <- normal(mu = -1.1, sigma = 2.1)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(2, 0))
  expect_equal(num$hyper, c(-1.1, 2.1, 0))
  expect_equal(num$hyper_names, c("mu", "sigma"))
  p <- normal(mu = 1.1, sigma = 2.1, square = TRUE)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(2, 1))
})

test_that("normal prior does not allow non-positive std", {
  reason <- "<sigma> must be positive! found ="
  expect_error(normal(mu = 1, sigma = -2), reason)
  expect_error(normal(mu = 1, sigma = 0), reason)
})

test_that("student-t prior is parsed correctly", {
  p <- student_t(20)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(3, 0))
  expect_equal(num$hyper, c(20, 0, 0))
  expect_equal(num$hyper_names, c("nu"))
})

test_that("gamma prior is parsed correctly", {
  p <- gam(shape = 2, inv_scale = 5)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(4, 0))
  expect_equal(num$hyper, c(2, 5, 0))
  expect_equal(num$hyper_names, c("alpha", "beta"))
})

test_that("inverse gamma prior is parsed correctly", {
  p <- igam(shape = 3, scale = 6, square = TRUE)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(5, 1))
  expect_equal(num$hyper, c(3, 6, 0))
  expect_equal(num$hyper_names, c("alpha", "beta"))
})

test_that("log-normal prior is parsed correctly", {
  p <- log_normal(3, 2)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(6, 0))
  expect_equal(num$hyper, c(3, 2, 0))
  expect_equal(num$hyper_names, c("mu", "sigma"))
})

test_that("beta prior is parsed correctly", {
  p <- bet(0.1, 0.5)
  expect_equal(p$alpha, 0.1)
  expect_equal(p$beta, 0.5)
})

test_that("invalid prior name or hyperparams cannot be given", {
  p <- list(dist = "stupid", square = TRUE)
  reason <- "given value 'stupid' for argument <distribution_name> is invalid"
  expect_error(prior_to_num(p), reason)
})

test_that("default prior warning works", {
  msg <- warn_msg_default_prior("param_desc", "param_name", "model_desc")
  expect_equal(nchar(msg), 187)
})
