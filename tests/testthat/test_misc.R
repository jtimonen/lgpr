library(lgpr)

# -------------------------------------------------------------------------

context("Argument validation")

test_that("check_allowed throws error for invalid input", {
  expect_error(check_allowed(arg = c("moi", "hei"), allowed = c("moi", "hei")))
  expect_error(check_allowed(arg = c("moi"), allowed = c("moi")))
  expect_error(check_allowed(arg = "juu", allowed = c("moi", "joo")))
  expect_error(check_allowed(arg = "juu", allowed = c("juu", "juu")))
})

test_that("check_allowed works correctly for valid input", {
  idx <- check_allowed(arg = c("hei"), allowed = c("moi", "hei"))
  expect_equal(idx, 2)
})

test_that("check_type works correctly", {
  out1 <- check_type(3.4, "numeric")
  out2 <- check_type(3.4, c("numeric", "character"))
  expect_true(out1)
  expect_true(out2)
  reason <- "has invalid type 'numeric'"
  expect_error(check_type(3.4, "list"), reason)
})

test_that("check_length works correctly", {
  out1 <- check_length("moi", 1)
  out2 <- check_length(c(4, 3, 2), 3)
  expect_true(out1)
  expect_true(out2)
  reason <- "has length 3, but its length should be 5"
  expect_error(check_length(c(4, 3, 2), 5), reason)
})

test_that("check_lengths works correctly", {
  a1 <- c(1, 2)
  a2 <- c(3, 2, 1)
  a3 <- c(3, 2, 1)
  expect_true(check_lengths(a2, a3))
  reason <- "lengths of a1 and a3 must match! found = 2 and 3"
  expect_error(check_lengths(a1, a3), reason)
})

test_that("check_not_null works correctly", {
  a1 <- c(1, 2)
  a2 <- NULL
  expect_true(check_not_null(a1))
  reason <- "a2 should not be NULL"
  expect_error(check_not_null(a2), reason)
})

test_that("check_interval works correctly", {
  x <- 0.3
  reason <- "<x> must be on the interval"
  expect_error(check_interval(x, 1, 2), reason)
  expect_true(check_interval(1.53, 1, 2))
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
