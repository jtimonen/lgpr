library(lgpr)

context("Input parsing: lgp_formula")

test_that("lgp_formula returns correct object", {
  f <- y ~ gp(x) + categorical(u) + mask(z) * gp(x) + gp_ns(age) + zerosum(q)
  a <- lgp_formula(f)
  expect_equal(.class2(a), "lgpformula")
  expect_equal(length(a@terms@summands), 5)
  expect_equal(length(rhs_variables(a@terms)), 6)
  expect_true(validObject(a))
  expect_gt(nchar(as.character(a)), 10)
})

test_that("quotes can be used in gp(), gp_ns() etc", {
  f <- y ~ gp("x") + categorical("u") + mask("z") * gp("x") +
    gp_ns("age") + zerosum("q")
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgp_formula works with a single lgpexpr", {
  f <- y ~ gp(x)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("response variable cannot be a covariate", {
  f <- y ~ gp(x) + mask(y) * gp(t)
  msg <- paste0(
    "invalid class “lgpformula” object: the response variable ",
    "cannot be also a covariate"
  )
  expect_error(lgp_formula(f), msg)
})

test_that("lgp_formula works with a single lgpterm", {
  f <- y ~ gp(x) * mask(z)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("two lgpterms can be summed", {
  f <- y ~ gp(x) * mask(z) + gp(x) * mask(z)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpterm and lgpexpr can be summed", {
  f <- y ~ gp(x) * mask(z) + mask("z")
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpexpr and lgpterm can be summed", {
  f <- y ~ gp(x) + mask("z") * gp_ns(aa)
  c <- .class2(lgp_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgp_formula throws error when input is not a formula", {
  expect_error(lgp_formula("a + b "))
})

test_that("lgp_formula throws error if invalid function or covariate", {
  expect_error(lgp_formula(y ~ notafunction(x)))
  expect_error(lgp_formula(y ~ gp("")))
  expect_error(lgp_formula(y ~ gp(x) * gp(y) * gp(z)))
  expect_error(lgp_formula(y ~ x + a))
  expect_error(lgp_formula(y ~ gp(x)(y)))
})
