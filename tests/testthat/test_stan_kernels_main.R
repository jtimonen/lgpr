library(lgpr)

# -------------------------------------------------------------------------

context("Stan high-level functions")

# Model
m <- create_model(y ~ zs(id) * gp(age) + het(id) * gp_vm(dis_age) +
  categ(sex) + gp(age) + gp(blood),
data = testdata_001
)

# Input
input <- m@stan_input
x_cat <- input$x_cat
x_mask <- input$x_cont_mask
num_levels <- input$x_cat_num_levels
comp <- input$components
n <- dim(x_cat)[2]
K_const <- kernel_const_all(
  n, n, x_cat, x_cat, x_mask, x_mask,
  num_levels, comp
)

# Params
alpha <- c(1, 1, 1, 1, 1)
ell <- c(1, 1, 1, 1)
beta <- c(0.5, 1.0)
x <- input$x_cont
x_unnorm <- input$x_cont_unnorm
vm_params <- input$vm_params
ix <- input$idx_expand

# All kernels
teff_zero <- dollar(input, "teff_zero")
K <- kernel_all(
  n, n, K_const, comp, x, x, x_unnorm, x_unnorm,
  alpha, ell, 0.5, beta, NaN,
  vm_params, ix, ix, teff_zero
)

test_that("kernel_const_all works correctly", {
  expect_equal(length(K_const), 5)
  expect_equal(dim(K_const[[1]]), c(24, 24))
  S1 <- sum(K_const[[1]])
  expect_lt(abs(S1), 1e-6)
  expect_equal(sum(K_const[[4]]), 24 * 24)
  expect_equal(sum(K_const[[5]]), 24 * 24)
})

test_that("kernel_all works correctly", {
  expect_equal(length(K), 5)
  expect_equal(dim(K[[1]]), c(24, 24))
})

test_that("gp_posterior works correctly", {
  y <- testdata_001$y
  fp <- gp_posterior(K, K, K, y, 1e-6, 1.0)
  expect_equal(length(fp), 12)

  # test that component-wise means sum to total mean
  f_sum <- fp[[1]]
  for (j in 2:5) {
    f_sum <- f_sum + fp[[j]]
  }
  diff <- f_sum - fp[[6]]
  expect_lt(max(abs(diff)), 1e-6)
})


# -------------------------------------------------------------------------

# Model
m <- create_model(y ~ age + sex, testdata_001)

# Input
input <- m@stan_input
x_cat <- input$x_cat
x_mask <- input$x_cont_mask
num_levels <- input$x_cat_num_levels
comp <- input$components
n <- dim(x_cat)[2]
K_const <- kernel_const_all(
  n, n, x_cat, x_cat, x_mask, x_mask,
  num_levels, comp
)

# Params
alpha <- c(1, 1)
ell <- c(1)
x <- input$x_cont
x_unnorm <- input$x_cont_unnorm
vm_params <- input$vm_params
ix <- input$idx_expand

# All kernels
teff_zero <- dollar(input, "teff_zero")
K <- kernel_all(
  n, n, K_const, comp, x, x, x_unnorm, x_unnorm,
  alpha, ell, NaN, NaN, NaN,
  vm_params, ix, ix, teff_zero
)

test_that("kernel_const_all works when no nonstationary components", {
  expect_equal(length(K_const), 2)
})

test_that("kernel_all works when no nonstationary components", {
  expect_equal(length(K), 2)
})
