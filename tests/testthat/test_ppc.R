library(lgpr)

# -------------------------------------------------------------------------

context("Prior predictive checks")

test_that("prior checks can be performed with gaussian obs model", {
  my_prior <- list(alpha = gam(2, 5), ell = gam(3, 5), sigma = igam(15, 10))
  suppressWarnings({
    y_rng <- prior_predict(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      prior = my_prior,
      data = testdata_001,
      iter = 200,
      chains = 2,
      refresh = 0,
      cores = 2
    )
    expect_equal(dim(y_rng), c(200, 24))
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
    sample_f = FALSE
  ), reason)

  suppressWarnings({
    y_rng <- prior_predict(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "nb",
      data = dat,
      iter = 200,
      chains = 1,
      refresh = 0,
    )
    d <- max(abs(y_rng - round(y_rng)))
    expect_lt(d, 1e-6)
  })
})


test_that("prior checks can be performed with beta bin obs model", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  suppressWarnings({
    y_rng <- prior_predict(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "bb",
      data = dat,
      iter = 1000,
      chains = 1,
      refresh = 0,
      num_trials = 10
    )
    d <- max(abs(y_rng - round(y_rng)))
    expect_lt(d, 1e-6)
  })
})
