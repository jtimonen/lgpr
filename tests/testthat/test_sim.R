library(lgpr)

# -------------------------------------------------------------------------

context("sim: simulating data")


test_that("gaussian data can be simulated", {
  expect_equal(
    length(names(simulate_data(
      N = 16,
      t_data = seq(6, 36, by = 6),
      covariates = c(2, 2),
      lengthscales = c(6, 6, 6, 6),
      relevances = c(1, 1, 1, 0),
      names = c("sex", "location"),
      t_jitter = 0.5
    ))),
    8
  )
})

test_that("poisson data can be simulated", {
  expect_equal(
    length(names(simulate_data(
      N = 10,
      t_data = seq(6, 36, by = 6),
      covariates = c(2, 2),
      noise_type = "poisson"
    ))),
    8
  )
})

test_that("neg binomial data can be simulated", {
  expect_equal(
    length(names(simulate_data(
      N = 4,
      t_data = c(1, 2, 3),
      covariates = c(2, 2),
      noise_type = "nb"
    ))),
    8
  )
})



# -------------------------------------------------------------------------

context("sim: kernel functions")
library(lgpr)

test_that("base kernels work correctly", {
  expect_equal(
    sim_kernel_zerosum(c(1, 2), c(3, 2, 1), M = 3, alpha = 1),
    matrix(c(-0.5, -0.5, 1.0, -0.5, 1.0, -0.5),
      nrow = 2, ncol = 3, byrow = TRUE
    )
  )
  expect_equal(
    sim_kernel_bin(c(1, 2), c(3, 2, 1), pos_class = 2),
    matrix(c(0, 0, 0, 1, 0, 0), nrow = 2, ncol = 3, byrow = FALSE)
  )
  expect_equal(
    sim_kernel_se(-2, -2, ell = 20),
    matrix(1)
  )
  expect_equal(
    dim(sim_kernel_ns(c(1, 1, 2), c(0, 1), ell = 1, a = 1, b = -10, c = 1)),
    c(3, 2)
  )
})

test_that("base kernels give errors when supposed to", {
  expect_error(sim_kernel_se(0, c(-1, 0, 1), ell = 0))
  expect_error(sim_kernel_cat(0, c(-1, 0, 1), ell = -3, alpha = 1))
  expect_error(sim_kernel_ns(0, c(-1, 0, 1),
    ell = 1, alpha = -1, a = 1, b = 0, c = 1
  ))
  expect_error(sim_kernel_ns(0, c(-1, 0, 1),
    alpha = 1, a = 1, b = 0, c = 1
  ))
})
