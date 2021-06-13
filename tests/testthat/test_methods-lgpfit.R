library(lgpr)
source("helpers/SW.R")

# -------------------------------------------------------------------------
SW({
  fit <- example_fit(refresh = 0, quiet = TRUE)
})

context("Methods for lgpfit objects")

test_that("postprocessing information can be removed and computed again", {
  expect_error(postproc(fit), "already contains postprocessing information")
  fit <- clear_postproc(fit)
  expect_true(!contains_postproc(fit))
  fit <- postproc(fit, verbose = FALSE)
  expect_true(contains_postproc(fit))
})
