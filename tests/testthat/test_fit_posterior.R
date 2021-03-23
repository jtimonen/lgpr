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
  expect_error(get_y_rng(fit), "y rng was not done!")
  expect_error(plot_warp(fit), "the model does not have warping parameters")
  expect_error(plot_beta(fit), "there are no heterogeneous effects")
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
  })
  expect_s4_class(fit, "lgpfit")
  expect_output(show(fit@model))
  r <- relevances(fit)
  expect_equal(length(r), 3)
  s <- select(fit)
  expect_equal(length(s$Component), 3)
  t <- seq(0, 40, by = 1)
  x_pred <- new_x(dat, t)
  p <- pred(fit, x_pred, reduce = NULL, draws = c(1:3), verbose = FALSE)
  plt1 <- plot_pred(fit, x_pred, p) # [0,1] scale
  plt2 <- plot_f(fit, x_pred, p)
  expect_s3_class(plt1, "ggplot")
  expect_s3_class(plt2, "ggplot")
})

test_that("verbose mode can be used", {
  dat <- testdata_001
  expect_output(
    suppressWarnings({
      lgp(
        formula = y ~ gp(age) + zs(sex) * gp(age),
        data = dat,
        iter = 20,
        chains = 1,
      )
    })
  )
})

test_that("everything works when the dependent variable is not named 'y'", {
  dat <- testdata_001
  colnames(dat)[6] <- "depvar"
  suppressWarnings({
    fit <- lgp(depvar ~ age + sex,
      data = dat, iter = 100, chains = 1,
      refresh = 0
    )
  })
  t <- seq(1, 50, by = 2)
  x_pred <- new_x(dat, t, x = "age")
  p <- pred(fit, x_pred, reduce = mean, verbose = FALSE)
  expect_equal(get_y_name(fit), "depvar")
  h1 <- plot_components(fit, x = x_pred, pred = p)
  h2 <- plot_f(fit, x = x_pred, pred = p)
  h3 <- plot_pred(fit, x = x_pred, pred = p, t_name = "age")
  expect_equal(length(h1), 3)
  expect_s3_class(h2, "ggplot")
  expect_s3_class(h3, "ggplot")
})


test_that("everything works when id variable is not named 'id'", {
  dat <- testdata_001
  dat$SUBJECT <- dat$id
  dat$id <- NULL
  suppressWarnings({
    fit <- lgp(y ~ gp(age) + gp_vm(dis_age),
      data = dat, iter = 100, chains = 1,
      refresh = 0
    )
  })

  # get_teff_obs works correctly
  expect_error(get_teff_obs(dat, x_ns = "dis_age"), "'id' not found in <data>")
  teff <- get_teff_obs(dat, x_ns = "dis_age", group_by = "SUBJECT")
  expect_equal(as.numeric(teff), c(21, 21, NaN, NaN))

  # need to give group_by = "SUBJECT" to every function now
  t <- seq(1, 50, by = 2)
  x_pred <- new_x(dat, t, x = "age", group_by = "SUBJECT", x_ns = "dis_age")
  p <- pred(fit, x_pred, reduce = mean, verbose = FALSE)
  h1 <- plot_components(fit, x = x_pred, pred = p, group_by = "SUBJECT")
  h2 <- plot_f(fit, x = x_pred, pred = p, group_by = "SUBJECT")
  h3 <- plot_pred(fit,
    x = x_pred, pred = p, t_name = "age",
    group_by = "SUBJECT"
  )
  expect_equal(length(h1), 3)
  expect_s3_class(h2, "ggplot")
  expect_s3_class(h3, "ggplot")
})
