library(lgpr)
source("helpers/SW.R")


# -------------------------------------------------------------------------

context("Computing predictions (f marginalized)")

N_ITER <- 20
N_CHAINS <- 1
DAT <- testdata_001


test_that("predictions can be computed (f marginalized)", {
  my_prior <- list(effect_time_info = list(
    lower = 20,
    upper = 30,
    backwards = FALSE,
    zero = 0
  ))

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

  # Create test points
  # x_pred <- new_x(data = DAT, x_values = seq(0, 40, 0.5))
  x_pred <- NULL
  p <- pred(fit, x_pred, verbose = FALSE)
  expect_equal(dim(p), c(6, 24)) #

  # TODO: make pred return (scaled!) y_pred and test it!!
})
