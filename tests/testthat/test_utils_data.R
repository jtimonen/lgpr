library(lgpr)

# -------------------------------------------------------------------------

context("Data utility functions")

test_that("computing observed effect times from data works correctly", {
  age <- c(10, 20, 30, 10, 20, 30)
  dage <- c(-10, 0, 10, NaN, NaN, NaN)
  id <- as.factor(c(4, 4, 4, 6, 6, 6))
  df <- data.frame(id, age, dage)
  t0 <- get_teff_obs(df, x_ns = "dage")
  expect_equal(names(t0), c("4", "6"))
  expect_equal(as.numeric(t0), c(20, NaN))
})

test_that("data validation works correctly", {
  dat <- testdata_001
  dat$sex[1] <- "Female"
  expect_error(validate_data(dat), "do not all have the same")
  expect_true(validate_data(testdata_001))
})

test_that("new_data works correctly without x_ns", {
  dat <- testdata_001
  x_new <- new_data(dat, x_values = seq(-1, 1, by = 0.5))
  expect_equal(dim(x_new), c(20, 3))
  expect_equal(colnames(x_new), c("id", "age", "sex"))
})

test_that("new_data works correctly with x_ns", {
  dat <- testdata_001
  x_new <- new_data(dat, x_values = seq(-1, 1, by = 0.5), x_ns = "dis_age")
  expect_equal(dim(x_new), c(20, 4))
  expect_equal(colnames(x_new), c("id", "age", "dis_age", "sex"))
  num_nans <- sum(is.nan(x_new$dis_age))
  expect_equal(num_nans, 10)
})
