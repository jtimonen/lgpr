library(lgpr)
tryCatch(
  {
    source("helpers/SW.R")
  },
  error = function(e) {
    cat("unable to source test helpers\n")
  }
)


# -------------------------------------------------------------------------

context("Computing predictions (f marginalized)")

N_ITER <- 20
N_CHAINS <- 1
DAT <- testdata_001

my_prior <- list(effect_time_info = list(
  lower = 20,
  upper = 30,
  backwards = FALSE,
  zero = 0
))

test_that("predictions can be computed (f marginalized)", {


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
    pred(fit, x_pred, verbose = FALSE),
    "variable 'dis_age' not found in <x>!"
  )

  # Compute predictions with proper x
  x_pred <- new_x(data = DAT, x_values = seq(0, 40, 0.5), x_ns = "dis_age")
  p <- pred(fit, x_pred, verbose = FALSE, reduce = NULL)
  expect_s4_class(p, "GaussianPrediction")

  # Compute predictions with same x as in data
  p1 <- pred(fit, DAT, verbose = FALSE)
  p2 <- pred(fit, NULL, verbose = FALSE)
  expect_equal(p1, p2)

  # TODO: make pred return (scaled!) y_pred and test it!!
})


context("Computing predictions (f latent)")

test_that("predictions can be computed (f latent)", {
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
