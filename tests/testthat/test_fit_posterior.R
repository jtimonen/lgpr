library(lgpr)

# -------------------------------------------------------------------------

context("Posterior sampling and optimization")

my_prior <- list(
  alpha = normal(1, 0.1),
  ell = normal(1, 0.1),
  sigma = normal(1, 0.1)
)

dat <- testdata_001
model <- create_model(y ~ gp(age), dat, prior = my_prior)

test_that("sample_model can do posterior sampling", {
  fit <- sample_model(
    model = model,
    iter = 1200,
    chains = 1,
    refresh = 0,
    seed = 123
  )

  expect_s4_class(fit, "lgpfit")
  p1 <- plot_draws(fit)
  p2 <- plot_draws(fit, type = "trace")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  r <- relevances(fit)
  expect_equal(sum(r$Relevance), 1.0)
})

test_that("optimize_model can optimize MAP parameters", {
  fit <- optimize_model(model, iter = 10, seed = 123)
  expect_equal(class(fit), "list")
})


# -------------------------------------------------------------------------

context("Posterior sampling of f with different obs models")

test_that("f can be sampled with gaussian likelihood", {
  dat <- testdata_001
  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex),
      sample_f = TRUE,
      data = dat,
      iter = 200,
      chains = 1,
      refresh = 0,
    )
    expect_s4_class(fit, "lgpfit")
    expect_equal(get_stan_input(fit)$is_f_sampled, 1)
  })
})

test_that("f can be sampled with poisson likelihood", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "poisson",
      data = dat,
      iter = 200,
      chains = 1,
      refresh = 0,
    )
    expect_s4_class(fit, "lgpfit")
  })
})

test_that("f can be sampled with nb likelihood", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + categ(sex) * gp(age),
      likelihood = "nb",
      data = dat,
      iter = 200,
      chains = 1,
      refresh = 0,
    )
    expect_s4_class(fit, "lgpfit")
    p1 <- plot_draws(fit)
    p2 <- plot_draws(fit, regex_pars = "f_latent")
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")

    a <- get_f(fit)
    expect_equal(dim(a$f$`gp(age)`), c(100, 24))
    r <- relevances(fit)
    expect_equal(dim(r), c(3, 2))
  })
})

test_that("f can be sampled with binomial likelihood", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + zs(sex) * gp(age),
      likelihood = "binomial",
      data = dat,
      iter = 200,
      chains = 1,
      refresh = 0,
      num_trials = 10
    )
    expect_s4_class(fit, "lgpfit")
  })
})

test_that("f can be sampled with beta-binomial likelihood", {
  dat <- testdata_001
  dat$y <- round(exp(dat$y))
  suppressWarnings({
    fit <- lgp(
      formula = y ~ gp(age) + zs(sex) * gp(age),
      likelihood = "bb",
      data = dat,
      iter = 200,
      chains = 1,
      refresh = 0,
      num_trials = 10
    )
    expect_s4_class(fit, "lgpfit")
    expect_error(gp_posteriors(fit, dat), "observation model must be gaussian")
  })
})
