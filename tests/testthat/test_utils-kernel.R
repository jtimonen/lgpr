set.seed(9393)
library(lgpr)
# source("tests/testthat/helpers/SW.R")
source("helpers/SW.R")
N_ITER <- 19
N_CHAINS <- 2

# -------------------------------------------------------------------------

context("Kernel matrix computations")

test_that("matrices have correct shape when model has multiple components", {
  J <- 4
  N_OBS <- 30
  SW({
    fit <- example_fit(iter = N_ITER, chains = N_CHAINS)
  })
  S <- fit@num_draws

  # No new x input and no reduce
  mats <- posterior_f(fit, reduce = NULL, debug_km = TRUE, verbose = FALSE)
  expect_equal(dim(mats$K), c(S, J, N_OBS, N_OBS))
  expect_equal(dim(mats$Ks), c(S, J, N_OBS, N_OBS))
  expect_equal(dim(mats$Kss), c(S, J, N_OBS, N_OBS))

  # With new x input and no reduce
  x_out <- new_x(fit@model@data, x_values = seq(0, 100, by = 2.5))
  N_OUT <- nrow(x_out)
  mats <- posterior_f(fit,
    x = x_out, reduce = NULL, draws = c(1, 3, 5),
    debug_km = TRUE, verbose = FALSE
  )
  expect_equal(dim(mats$K), c(3, J, N_OBS, N_OBS))
  expect_equal(dim(mats$Ks), c(3, J, N_OUT, N_OBS))
  expect_equal(dim(mats$Kss), c(3, J, N_OUT, N_OUT))

  # With new x input and reduce (mean params)
  mats <- posterior_f(fit, x = x_out, debug_km = TRUE, verbose = FALSE)
  expect_equal(dim(mats$K), c(1, J, N_OBS, N_OBS))
  expect_equal(dim(mats$Ks), c(1, J, N_OUT, N_OBS))
  expect_equal(dim(mats$Kss), c(1, J, N_OUT, N_OUT))

  # Expect output in debug mode
  expect_output(
    posterior_f(fit, debug_km = TRUE, debug_dims = TRUE, verbose = FALSE)
  )

  # Expect output in verbose mode
  expect_output(
    posterior_f(fit, debug_km = TRUE, debug_dims = FALSE, verbose = TRUE)
  )
})


test_that("matrices have correct shape when model has only one component", {
  J <- 1
  N_OBS <- 30
  SW({
    fit <- example_fit(formula = y ~ id, iter = N_ITER, chains = N_CHAINS)
  })
  S <- fit@num_draws

  # No new x input and no reduce
  mats <- posterior_f(fit, reduce = NULL, debug_km = TRUE, verbose = FALSE)
  expect_equal(dim(mats$K), c(S, J, N_OBS, N_OBS))
  expect_equal(dim(mats$Ks), c(S, J, N_OBS, N_OBS))
  expect_equal(dim(mats$Kss), c(S, J, N_OBS, N_OBS))

  # With new x (only one point) input and no reduce
  x_out <- new_x(fit@model@data, x_values = 3.2)
  N_OUT <- nrow(x_out)
  mats <- posterior_f(fit,
    x = x_out, reduce = NULL, draws = c(1, 3, 5),
    debug_km = TRUE, verbose = FALSE
  )
  expect_equal(dim(mats$K), c(3, J, N_OBS, N_OBS))
  expect_equal(dim(mats$Ks), c(3, J, N_OUT, N_OBS))
  expect_equal(dim(mats$Kss), c(3, J, N_OUT, N_OUT))

  # With new x input and reduce
  x_out <- new_x(fit@model@data, x_values = c(4.5, 7.7, 90.9))
  N_OUT <- nrow(x_out)
  mats <- posterior_f(fit, x = x_out, debug_km = TRUE, verbose = FALSE)
  expect_equal(dim(mats$K), c(1, J, N_OBS, N_OBS))
  expect_equal(dim(mats$Ks), c(1, J, N_OUT, N_OBS))
  expect_equal(dim(mats$Kss), c(1, J, N_OUT, N_OUT))
})
