context("kernel functions")
library(lgpr)

test_that("base kernels work correctly", {
  expect_equal(
    kernel_zerosum(c(1,2),c(3,2,1), M = 3, alpha = 1),
    matrix(c(-0.5,-0.5,1.0,-0.5,1.0,-0.5), nrow = 2, ncol = 3, byrow = TRUE)
  )
  expect_equal(
    kernel_bin(c(1,2),c(3,2,1), pos_class = 2),
    matrix(c(0,0,0,1,0,0), nrow = 2, ncol = 3, byrow = FALSE)
  )
  expect_equal(
    kernel_se(-2,-2,ell=20),
    matrix(1)
  )
  expect_equal(
    dim(kernel_ns(c(1,1,2),c(0,1),ell=1,a=1,b=-10,c=1)),
    c(3,2)
  )
})

test_that("base kernels give errors when supposed to", {
  expect_error(kernel_se(0,c(-1,0,1),ell=0))
  expect_error(kernel_cat(0,c(-1,0,1),ell=-3, alpha=1))
  expect_error(kernel_ns(0,c(-1,0,1), ell=1, alpha=-1, a = 1, b = 0, c = 1))
  expect_error(kernel_ns(0,c(-1,0,1), alpha=1, a = 1, b = 0, c = 1))
})
