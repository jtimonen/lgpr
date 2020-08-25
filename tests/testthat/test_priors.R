library(lgpr)

# -------------------------------------------------------------------------

context("Prior parsing")


test_that("uniform prior is parsed correctly", {
  suppressWarnings({
    p <- uniform()
    num <- prior_to_num(p)
    expect_equal(num$prior, c(1, 0))
  })
})

test_that("defining uniform prior gives warning", {
  expect_warning(uniform())
})

test_that("normal prior is parsed correctly", {
  p <- normal(mu = -1.1, sigma = 2.1)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(2, 0))
  expect_equal(num$hyper, c(-1.1, 2.1, 0))
  expect_equal(num$hyper_names, c("mu", "sigma"))
  p <- normal(mu = 1.1, sigma = 2.1, square = TRUE)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(2, 1))
})

test_that("normal prior does not allow non-positive std", {
  reason <- "<sigma> must be positive! found ="
  expect_error(normal(mu = 1, sigma = -2), reason)
  expect_error(normal(mu = 1, sigma = 0), reason)
})

test_that("student-t prior is parsed correctly", {
  p <- student_t(20)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(3, 0))
  expect_equal(num$hyper, c(20, 0, 0))
  expect_equal(num$hyper_names, c("nu"))
})

test_that("gamma prior is parsed correctly", {
  p <- gam(shape = 2, inv_scale = 5)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(4, 0))
  expect_equal(num$hyper, c(2, 5, 0))
  expect_equal(num$hyper_names, c("alpha", "beta"))
})

test_that("inverse gamma prior is parsed correctly", {
  p <- igam(shape = 3, scale = 6, square = TRUE)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(5, 1))
  expect_equal(num$hyper, c(3, 6, 0))
  expect_equal(num$hyper_names, c("alpha", "beta"))
})

test_that("log-normal prior is parsed correctly", {
  p <- log_normal(3, 2)
  num <- prior_to_num(p)
  expect_equal(num$prior, c(6, 0))
  expect_equal(num$hyper, c(3, 2, 0))
  expect_equal(num$hyper_names, c("mu", "sigma"))
})

test_that("beta prior is parsed correctly", {
  p <- bet(0.1, 0.5)
  expect_equal(p$alpha, 0.1)
  expect_equal(p$beta, 0.5)
})

test_that("invalid prior name or hyperparams cannot be given", {
  p <- list(dist = "stupid", square = TRUE)
  reason <- "given value 'stupid' for argument <distribution_name> is invalid"
  expect_error(prior_to_num(p), reason)
})
