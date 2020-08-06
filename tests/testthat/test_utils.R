library(lgpr)

# -------------------------------------------------------------------------

context("Argument checker")

test_that("argument checker throws error for invalid input", {
  expect_error(argument_check(arg = c("moi", "hei"), allowed=c("moi", "hei") ))
  expect_error(argument_check(arg = c("moi"), allowed=c("moi") ))
  expect_error(argument_check(arg = "juu", allowed=c("moi", "joo") ))
  expect_error(argument_check(arg = "juu", allowed=c("juu", "juu") ))
})

test_that("argument checker works correctly for valid input", {
  idx <- argument_check(arg = c("hei"), allowed=c("moi", "hei") )
  expect_equal(idx, 2)
})

# -------------------------------------------------------------------------

context("Inverse-gamma distribution")

test_that("inverse gamma density works", {
  a <- 2
  b <- 3
  diff <- abs(dinvgamma_stanlike(1, a, b) - 0.4480836)
  expect_lt(diff, 1e-6)
})

test_that("inverse gamma quantile works", {
  a <- 2
  b <- 3
  diff <- abs(qinvgamma_stanlike(0.9, a, b) - 5.641095)
  expect_lt(diff, 1e-6)
})

test_that("plotting the inverse gamma distribution works", {
  a <- 2
  b <- 3
  plot <- plot_invgamma(a, b)
  class_str <- as.character(class(plot))
  expect_equal(class_str, c('gg', 'ggplot'))
})

# -------------------------------------------------------------------------

context("Misc utility functions")

test_that("repvec works correctly", {
  expect_equal(repvec(3, 1), matrix(3))
  expect_equal(repvec(3, 3), matrix(c(3, 3, 3), ncol = 1, nrow = 3))
  expect_equal(
    repvec(c(2, 3, 4), 2),
    matrix(c(2, 3, 4, 2, 3, 4), ncol = 3, nrow = 2, byrow = TRUE)
  )
})

test_that("get_stan_model returns a stanmodel", {
  sm <- get_stan_model()
  expect_equal(as.character(class(sm)), 'stanmodel')
})

test_that("example_call returns text", {
  call_str <- example_call()
  expect_equal(as.character(class(call_str)), 'character')
  expect_gt(nchar(call_str), 100)
})

test_that("color_set works", {
  col1 <- color_set('blue')
  col2 <- color_set('red')
  expect_equal(nchar(col1), 7)
  expect_equal(nchar(col2), 7)
  expect_error(color_set('definitely_not_a_color_name'))
})

# -------------------------------------------------------------------------
