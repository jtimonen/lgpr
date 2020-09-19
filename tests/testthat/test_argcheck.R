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

test_that("check_all_leq works correctly", {
  x <- c(4, 3, 4, 5)
  y <- c(4, 3, 2, 1)
  reason <- "value of <x> is larger than value of <y> at index 3"
  expect_error(check_all_leq(x, y), reason)
  expect_true(check_all_leq(c(1, 1, 1, 1), y))
})

test_that("check_not_named works correctly", {
  x <- c(4, 3, 4, 5)
  expect_true(check_not_named(x))
  names(x) <- c("hei", "hey", "ho", "jea")
  reason <- "<x> should not have names! found = "
  expect_error(check_not_named(x), reason)
})

test_that("check_named works correctly", {
  x <- c(4, 3, 4, 5)
  expect_error(check_named(x), "<x> must have names")
  names(x) <- c("hei", "hey", "ho", "jea")
  expect_true(check_named(x))
})


test_that("check_numeric works correctly", {
  reason <- "must be numeric"
  expect_error(check_numeric("moi"), reason)
  expect_true(check_numeric(1))
})

test_that("check_null works correctly", {
  a <- 123
  msg <- "should be NULL! Reason: no reason"
  expect_error(check_null(a, "no reason"), msg)
  expect_true(check_null(NULL))
})

test_that("check_false works correctly", {
  a <- TRUE
  reason <- "to be FALSE"
  expect_error(check_false(a), reason)
  expect_error(check_false(2), reason)
  expect_true(check_false(0))
  expect_true(check_false(!a))
})

test_that("check_dim works correctly", {
  a <- array(0, c(2, 3, 4))
  expect_true(check_dim(a, 3))
  reason <- "number of dimensions of <a> must be 2! found = 3"
  expect_error(check_dim(a, 2), reason)
})
