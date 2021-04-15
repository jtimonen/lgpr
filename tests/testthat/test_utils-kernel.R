set.seed(93391)
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
  expect_equal(length(mats$K), S)
  tmp <- mats$Ks[[1]]
  expect_equal(dim(tmp[[1]]), c(N_OBS, N_OBS))
  expect_equal(length(mats$Kss[[1]]), J)

  # With new x input and no reduce
  x_out <- new_x(fit@model@data, x_values = seq(0, 100, by = 2.5))
  N_OUT <- nrow(x_out)
  mats <- posterior_f(fit,
    x = x_out, reduce = NULL, draws = c(1, 3, 5),
    debug_km = TRUE, verbose = FALSE
  )
  expect_equal(length(mats$K), 3)
  tmp <- mats$Ks[[1]]
  expect_equal(dim(tmp[[1]]), c(N_OUT, N_OBS))
  expect_equal(length(mats$Kss[[1]]), J)

  # With new x input and reduce (mean params)
  mats <- posterior_f(fit, x = x_out, debug_km = TRUE, verbose = FALSE)
  expect_equal(length(mats$K), 1)
  tmp <- mats$Ks[[1]]
  expect_equal(dim(tmp[[1]]), c(N_OUT, N_OBS))
  expect_equal(length(mats$Kss[[1]]), J)

  # Expect output in debug mode
  expect_message(
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
  expect_equal(length(mats$K), S)
  tmp <- mats$Ks[[1]]
  expect_equal(dim(tmp[[1]]), c(N_OBS, N_OBS))
  expect_equal(length(mats$Kss[[1]]), J)

  # With new x (only one point) input and no reduce
  x_out <- new_x(fit@model@data, x_values = 3.2)
  N_OUT <- nrow(x_out)
  mats <- posterior_f(fit,
    x = x_out, reduce = NULL, draws = c(1, 3, 5),
    debug_km = TRUE, verbose = FALSE
  )
  expect_equal(length(mats$K), 3)
  s_idx <- 2
  j_idx <- sample.int(n = J, size = 1)
  tmp <- mats$Ks[[s_idx]]
  expect_equal(dim(tmp[[j_idx]]), c(N_OUT, N_OBS))
  expect_equal(length(mats$Kss[[s_idx]]), J)

  # With new x input and reduce
  x_out <- new_x(fit@model@data, x_values = c(4.5, 7.7, 90.9))
  N_OUT <- nrow(x_out)
  mats <- posterior_f(fit, x = x_out, debug_km = TRUE, verbose = FALSE)
  expect_equal(length(mats$K), 1)
  s_idx <- 1
  j_idx <- sample.int(n = J, size = 1)
  tmp <- mats$Ks[[s_idx]]
  expect_equal(dim(tmp[[j_idx]]), c(N_OUT, N_OBS))
  expect_equal(length(mats$Kss[[s_idx]]), J)
})
