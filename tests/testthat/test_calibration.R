library(lgpr)

# -------------------------------------------------------------------------

context("Simulation-based calibration")

test_that("sbc can be performed with gaussian obs model", {
  my_prior <- list(alpha = gam(2, 5), ell = gam(3, 5), sigma = igam(15, 10))
  suppressWarnings({
    fit <- calibrate(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      prior = my_prior,
      data = testdata_001,
      iter = 600,
      chains = 2,
      refresh = 0,
      cores = 2
    )
    expect_s4_class(fit, "lgpfit")
    p <- plot_draws(fit)
    expect_s3_class(p, "ggplot")
  })
})

test_that("sbc can be performed with neg bin obs model", {
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
    fit <- calibrate(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "nb",
      data = dat,
      iter = 600,
      chains = 1,
      refresh = 0,
    )
    expect_s4_class(fit, "lgpfit")
    p1 <- plot_draws(fit)
    p2 <- plot_draws(fit, regex_pars = "f_latent")
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
  })
})
