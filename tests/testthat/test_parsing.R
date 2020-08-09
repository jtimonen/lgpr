library(lgpr)


# -------------------------------------------------------------------------

context("Input parsing: parsing formulas")

test_that("parse_formula returns correct object", {
  f <- y ~ gp(x) + categorical(u) + mask(z) * gp(x) + gp_ns(age) + zerosum(q)
  a <- parse_formula(f)
  expect_equal(.class2(a), "lgpformula")
  expect_equal(length(a@terms@summands), 5)
  expect_equal(length(rhs_variables(a@terms)), 6)
  expect_true(validObject(a))
  expect_gt(nchar(as.character(a)), 10)
})

test_that("quotes can be used in gp(), gp_ns() etc", {
  f <- y ~ gp("x") + categorical("u") + mask("z") * gp("x") +
    gp_ns("age") + zerosum("q")
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("parse_formula works with a single lgpexpr", {
  f <- y ~ gp(x)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("response variable cannot be a covariate", {
  f <- y ~ gp(x) + mask(y) * gp(t)
  msg <- "the response variable cannot be also a covariate"
  expect_error(parse_formula(f), msg)
})

test_that("parse_formula works with a single lgpterm", {
  f <- y ~ gp(x) * mask(z)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("two lgpterms can be summed", {
  f <- y ~ gp(x) * mask(z) + gp(x) * mask(z)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpterm and lgpexpr can be summed", {
  f <- y ~ gp(x) * mask(z) + mask("z")
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("lgpexpr and lgpterm can be summed", {
  f <- y ~ gp(x) + mask("z") * gp_ns(aa)
  c <- .class2(parse_formula(f))
  expect_equal(c, "lgpformula")
})

test_that("parse_formula throws error when input is not a formula", {
  expect_error(parse_formula("a + b "))
})

test_that("parse_formula throws error if invalid function or covariate", {
  expect_error(parse_formula(y ~ notafunction(x)))
  expect_error(parse_formula(y ~ gp("")))
  expect_error(parse_formula(y ~ gp(x) * gp(y) * gp(z)))
  expect_error(parse_formula(y ~ x + a))
  expect_error(parse_formula(y ~ gp(x)(y)))
})

# -------------------------------------------------------------------------

context("Input parsing: parsing options")

test_that("parse_options does not need arguments", {
  a <- parse_options()
  expect_equal(names(a), c("verbose","skip_generated","sample_f","delta"))
})

test_that("parse_disease_options does not need arguments", {
  a <- parse_disease_options()
  expect_equal(names(a), c("uncertain","heterogeneous","vm_params"))
})

test_that("vm_params gets correct dimensions", {
  a <- parse_disease_options()
  b <- parse_disease_options(list(uncertain=TRUE, vm_params=NA))
  expect_equal(dim(a$vm_params), c(1,2))
  expect_equal(dim(b$vm_params), c(0,2))
})



# -------------------------------------------------------------------------

context("Input parsing: parsing likelihood")

test_that("parse_likelihood can be used", {
  a <- parse_likelihood('gaussian', NULL, NULL, c(1,2))
  expect_equal(names(a), c("likelihood","num_trials","c_hat"))
})

test_that("parse_likelihood errors correctly", {
  expect_error(parse_likelihood('gaussian', 0, NULL, c(1,2)))
  expect_error(parse_likelihood('gaussian', NULL, 0, c(1,2)))
})

test_that("parse_likelihood works correctly with binomial likelihood", {
  a <- parse_likelihood('binomial', 1.2, 9, c(1,2))
  b <- parse_likelihood('binomial', NULL, 9, c(1,2))
  expect_equal(a$likelihood, 4)
  expect_equal(length(a$num_trials), 2)
  diff <- abs(b$c_hat - (-1.609438))
  expect_lt(max(diff), 1e-6)
  expect_error(parse_likelihood('binomial', NULL, c(1,1,1), c(1,2)))
})
