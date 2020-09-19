library(lgpr)
set.seed(123)

# -------------------------------------------------------------------------

context("Data splitting functions")
dat <- testdata_001

test_that("data can be split by id", {
  s <- split_by_factor(dat, var_name = "id", test = 1)
  e <- c("train", "test", "i_train", "i_test")
  expect_equal(names(s), e)
  expect_equal(length(s$i_test), 6)
  expect_equal(length(s$i_train), 24 - 6)
})

test_that("data can be split by time point", {
  s <- split_within_factor(dat, var_name = "id", idx_test = 3)
  expect_equal(length(s$i_test), 4)
  expect_equal(length(s$i_train), 24 - 4)
})

test_that("data can be split by time point randomly", {
  s <- split_within_factor_random(dat, var_name = "id", k_test = 2)
  expect_equal(length(s$i_test), 8)
  expect_equal(length(s$i_train), 24 - 8)
})

test_that("data can be split totally randomly", {
  s <- split_random(dat, p_test = 0.5)
  expect_equal(length(s$i_test), 12)
  expect_equal(length(s$i_train), 24 - 12)
})


# -------------------------------------------------------------------------

context("Other user assist functions")
dat <- testdata_001

test_that("adding a factor works correctly", {
  country <- c("FIN", "FIN", "EST", "EST")
  expect_error(add_factor(dat, country), "<x> must have names")
  names(country) <- c(1, 3, 2, 4)
  newdat <- add_factor(dat, country)
  expect_true("country" %in% names(newdat))
  expect_error(add_factor(newdat, country), "already contains")
})

test_that("adding disease age works correctly", {
  t_init <- c(10, 10, 10, 10)
  names(t_init) <- c(1, 2, 3, 4)
  expect_error(add_dis_age(dat, t_init), "already contains")

  newdat <- dat
  newdat$dis_age <- NULL
  newdat <- add_dis_age(newdat, t_init)
  expect_true("dis_age" %in% names(newdat))
})

test_that("adjusted_c_hat works correctly", {
  reason <- "inputs must have same length"
  expect_error(adjusted_c_hat(c(1, 1), norm_factors = c(1, 2, 3)), reason)
  reason <- "must be all positive"
  expect_error(adjusted_c_hat(c(1, 1, 0), norm_factors = c(1, 2, 0)), reason)
  reason <- "must have only integer values"
  expect_error(adjusted_c_hat(c(1, 1.2, 0), norm_factors = c(1, 2, 3)), reason)
  c_hat <- adjusted_c_hat(c(1, 1, 1), norm_factors = c(1, 2, 3))
  diff <- abs(c_hat - c(0.0000000, 0.6931472, 1.0986123))
  expect_lt(max(diff), 1e-6)
})
