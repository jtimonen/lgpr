library(lgpr)
source("helpers/SW.R")

# -------------------------------------------------------------------------

context("Prior predictive checks")

test_that("prior checks can be performed with gaussian obs model", {
  my_prior <- list(alpha = gam(2, 5), ell = gam(3, 5), sigma = igam(15, 10))
  dat <- testdata_001
  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      prior = my_prior,
      data = dat,
      iter = 200,
      chains = 2,
      refresh = 0,
      cores = 2,
      prior_only = TRUE,
      sample_f = TRUE,
      options = list(do_yrng = TRUE)
    )
    p <- ppc(fit, testdata_001)
    expect_s3_class(p, "ggplot")
  })
})

test_that("prior checks can be performed with neg bin obs model", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  reason <- "sample_f must be TRUE when likelihood is nb"
  expect_error(create_model(
    formula = y ~ gp(age) + categ(sex) * gp(age),
    likelihood = "nb",
    data = dat,
    sample_f = FALSE,
  ), reason)

  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "nb",
      data = dat,
      iter = 200,
      chains = 1,
      refresh = 0,
      prior_only = TRUE,
      options = list(do_yrng = TRUE)
    )
    y_rng <- get_y_rng(fit)
    d <- max(abs(y_rng - round(y_rng)))
    expect_lt(d, 1e-6)
  })
})


test_that("prior checks can be performed with beta bin obs model", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "bb",
      data = dat,
      iter = 1000,
      chains = 1,
      refresh = 0,
      num_trials = 20,
      prior_only = TRUE,
      options = list(do_yrng = TRUE)
    )
    p <- ppc(fit, dat)
    expect_s3_class(p, "ggplot")
    y_rng <- get_y_rng(fit)
    d <- max(abs(y_rng - round(y_rng)))
    expect_lt(d, 1e-6)
  })
})
