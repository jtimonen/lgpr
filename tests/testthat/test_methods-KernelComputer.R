library(lgpr)
source("helpers/SW.R")

# -------------------------------------------------------------------------
dat <- testdata_001

context("Methods for KernelComputer objects")

test_that("a KernelComputer can be created and printed", {
  SW({
    f <- example_fit()
  })
  sf <- get_stanfit(f)
  x <- get_data(f)

  kc <- create_kernel_computer(
    f@model, sf, x, TRUE, NULL, NULL, get_stream()
  )
  expect_output(print(kc), "full_covariance = FALSE")
})
