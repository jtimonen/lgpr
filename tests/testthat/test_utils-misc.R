library(lgpr)

# -------------------------------------------------------------------------
context("Miscellaneous utility functions")

test_that("link_inv functions are inverses of the link functions", {
  x <- c(1.2, -1, 0, 2)
  for (likelihood in likelihood_list()) {
    y <- link_inv(x, likelihood)
    expect_equal(x, link(y, !!likelihood))
  }
})
