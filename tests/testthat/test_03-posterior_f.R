library(lgpr)
source("helpers/SW.R")
# source("tests/testthat/helpers/SW.R")

# -------------------------------------------------------------------------

context("Computing analytical posterior distribution of f")

N_ITER <- 20
N_CHAINS <- 1
DAT <- testdata_001

my_prior <- list(effect_time_info = list(
  lower = 20,
  upper = 30,
  backwards = FALSE,
  zero = 0
))

test_that("posterior_f works (f marginalized)", {


  # Fit a model
  SW({
    fit <- lgp(
      y ~ gp(age) + het(id) * unc(id) * gp_vm(dis_age),
      data = DAT,
      iter = N_ITER,
      chains = N_CHAINS,
      refresh = 0,
      prior = my_prior
    )
  })

  # Try with misformatted x
  x_pred <- new_x(data = DAT, x_values = seq(0, 40, 0.5))
  expect_error(
    posterior_f(fit, x_pred, verbose = FALSE),
    "variable 'dis_age' not found in <x>!"
  )

  # Compute conditional function posteriors with proper x
  x_pred <- new_x(data = DAT, x_values = seq(0, 40, 0.5), x_ns = "dis_age")
  fp <- posterior_f(fit, x_pred, verbose = FALSE, reduce = NULL)
  expect_s4_class(fp, "FunctionPosteriors")
  expect_output(show(fp))

  # Plotting
  p1 <- plot(fp, verbose = FALSE)
  p2 <- plot(fp, group_by = "id", color_by = "sex", verbose = FALSE)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")

  # Same thing but with reduce to mean and x_pred = NULL
  fp <- posterior_f(fit, verbose = FALSE)
  expect_s4_class(fp, "FunctionPosteriors")
  p3 <- plot(fp, color_by = "sex")
  p4 <- plot(fp, group_by = "id", color_by = "dis_age")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})


# -------------------------------------------------------------------------

context("Getting draws from posterior of f")

test_that("posterior_f works (f latent)", {
  NEWDAT <- DAT
  NEWDAT$y <- round(exp(NEWDAT$y))

  # Fit a model
  SW({
    fit <- lgp(
      y ~ zs(id) + het(id) * unc(id) * gp_vm(dis_age),
      data = NEWDAT,
      iter = N_ITER,
      chains = N_CHAINS,
      refresh = 0,
      prior = my_prior,
      likelihood = "nb"
    )
  })

  # Compute predictions (TODO)
  x_pred <- new_x(data = DAT, x_values = seq(0, 60, 1.0), x_ns = "dis_age")
  # p <- pred(fit, x_pred, verbose = FALSE)
})
