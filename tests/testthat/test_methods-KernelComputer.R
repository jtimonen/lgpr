library(lgpr)
source("helpers/SW.R")

# -------------------------------------------------------------------------
dat <- testdata_001

context("Methods for KernelComputer objects")

test_that("a KernelComputer can be created and printed", {
  SW({
    f <- example_fit()
  })
  mod <- create_model(y ~ gp(age) + het(id) * gp_vm(dis_age), dat,
    prior = list(wrp = igam(14, 3))
  )
  kc <- create_kernel_computer(
    f@model, f@stan_fit,
    NULL, NULL, NULL, get_stream()
  )
  expect_output(print(kc), "Three same matrices: TRUE")
})
