library(lgpr)

# -------------------------------------------------------------------------

context("Argument checker")

test_that("argument checker throws error for invalid input", {
  expect_error(argument_check(arg = c("moi", "hei"), allowed = c("moi", "hei")))
  expect_error(argument_check(arg = c("moi"), allowed = c("moi")))
  expect_error(argument_check(arg = "juu", allowed = c("moi", "joo")))
  expect_error(argument_check(arg = "juu", allowed = c("juu", "juu")))
})

test_that("argument checker works correctly for valid input", {
  idx <- argument_check(arg = c("hei"), allowed = c("moi", "hei"))
  expect_equal(idx, 2)
})

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

context("Misc utility functions")

test_that("get_stan_model returns a stanmodel", {
  sm <- get_stan_model()
  expect_equal(as.character(class(sm)), "stanmodel")
})

test_that("color_set works", {
  col1 <- color_set("blue")
  col2 <- color_set("red")
  expect_equal(nchar(col1), 7)
  expect_equal(nchar(col2), 7)
  expect_error(color_set("definitely_not_a_color_name"))
})

# -------------------------------------------------------------------------
