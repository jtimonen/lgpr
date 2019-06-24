context("generating artificial data")
library(lgpr)

test_that("gaussian data can be simulated", {
  expect_equal(
    length(names(simulate_data(N           = 16,
                  t_data       = seq(6, 36, by = 6),
                  covariates   = c(    2,2),
                  lengthscales = c(6,6,6,6),
                  relevances   = c(1,1,1,0),
                  names        = c("sex", "location"),
                  t_jitter     = 0.5))),
   8
  )
})

test_that("poisson data can be simulated", {
  expect_equal(
    length(names(simulate_data(N           = 10,
                               t_data       = seq(6, 36, by = 6),
                               covariates   = c(    2,2),
                               noise_type   = "Poisson"))),
    8
  )
})

test_that("neg binomial data can be simulated", {
  expect_equal(
    length(names(simulate_data(N           = 4,
                               t_data       = c(1,2,3),
                               covariates   = c(    2,2),
                               noise_type   = "NB"))),
    8
  )
})
