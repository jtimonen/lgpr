library(lgpr)

# -------------------------------------------------------------------------

context("Inverse-gamma distribution")

test_that("inverse gamma density works and errors correctly", {
  a <- 2
  b <- 3
  diff <- abs(dinvgamma_stanlike(1, a, b) - 0.4480836)
  expect_lt(diff, 1e-6)
  expect_error(dinvgamma_stanlike(1, -1, 2))
  expect_error(dinvgamma_stanlike(1, 2, -1))
})

test_that("inverse gamma quantile works and errors correctly", {
  a <- 2
  b <- 3
  diff <- abs(qinvgamma_stanlike(0.9, a, b) - 5.641095)
  expect_lt(diff, 1e-6)
  expect_error(qinvgamma_stanlike(0.9, -1, 2))
  expect_error(qinvgamma_stanlike(0.9, 2, -1))
  expect_error(qinvgamma_stanlike(-0.1, 2, 2))
  expect_error(qinvgamma_stanlike(1.5, 2, 2))
})

test_that("plotting the inverse gamma distribution works", {
  a <- 2
  b <- 3
  p1 <- plot_invgamma(a, b)
  p2 <- plot_invgamma(a, b, return_quantiles = TRUE)
  c1 <- as.character(class(p1))
  c2 <- as.character(class(p2$plot))
  expect_equal(c1, c("gg", "ggplot"))
  expect_equal(c2, c("gg", "ggplot"))
  expect_error(plot_invgamma(a, b, IQR = 1.5))
})

# -------------------------------------------------------------------------

context("Utility functions")

test_that("ensure_len works correctly", {
  reason <- "length of <a> was expected to be 1 or 5, but found length 3"
  a <- c(1, 2, 3)
  expect_error(ensure_len(a, 5), reason)
  expect_equal(ensure_len(a, 3), a)
  expect_equal(ensure_len(4.2, 3), c(4.2, 4.2, 4.2))
})

test_that("dollar errors correctly", {
  a <- list(moi = c(2, 3, 4), juu = 1)
  reason <- paste0(
    "Element with name 'joo' not found in <a>! ",
    "Found elements:"
  )
  expect_error(dollar(a, "joo"), reason)
})

test_that("get_stan_model returns a stanmodel", {
  sm <- get_stan_model()
  expect_equal(as.character(class(sm)), "stanmodel")
})


test_that("add_sum_arraylist works correctly", {
  a <- 2 * diag(3)
  b <- 0.5 * diag(3)
  r <- 2 * diag(2)
  x1 <- list(a, b)
  x2 <- list(a, b, r)
  expect_error(add_sum_arraylist(x2), "non-conformable arrays")
  x <- add_sum_arraylist(x1)
  expect_equal(length(x), 3)
  expect_equal(x[[3]][1, 1], 2.5)
})
