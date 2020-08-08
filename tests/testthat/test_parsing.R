library(lgpr)

context("Input parsing: lgp_formula")

test_that("lgp_formula returns an object of class lgpformula", {
  f <- y ~ gp(x) + categorical(u) + mask(z) * gp(x) + gp_ns(age) + zerosum(q)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("quotes can be used in gp(), gp_ns() etc", {
  f <- y ~ gp("x") + categorical("u") + mask("z") * gp("x") +
    gp_ns("age") + zerosum("q")
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgp_formula works with a single lgpterm", {
  f <- y ~ gp(x)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgp_formula works with a single lgpproduct", {
  f <- y ~ gp(x) * mask(z)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("two lgpproducts can be summed", {
  f <- y ~ gp(x) * mask(z) + gp(x) * mask(z)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpproduct and lgpterm can be summed", {
  f <- y ~ gp(x) * mask(z) + mask("z")
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpterm and lgpproduct can be summed", {
  f <- y ~ gp(x) + mask("z") * gp_ns(aa)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgp_formula throws error when input is not a formula", {
  expect_error(lgp_formula("a + b "))
})
