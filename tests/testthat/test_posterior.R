library(lgpr)

# -------------------------------------------------------------------------

context("Function posteriors and posterior predictive distributions")

my_prior <- list(
  alpha = normal(1, 0.1),
  ell = normal(1, 0.1),
  sigma = normal(1, 0.1)
)

dat <- testdata_001
model <- create_model(y ~ gp(age) + zs(sex), dat, prior = my_prior)
fit <- sample_model(
  model = model,
  iter = 1200,
  chains = 1,
  refresh = 0,
  seed = 123
)
t <- seq(0, 50, by = 5)
X_pred <- new_x(dat, t)

test_that("new_x works correctly", {
  expect_equal(dim(X_pred), c(44, 3))
})

test_that("kernel parameter samples can be extracted", {
  theta <- get_draws_kernel_pars(fit)
  expect_equal(names(theta), c("alpha", "ell", "wrp", "beta", "teff"))
  expect_equal(nrow(theta$alpha), nrow(theta$teff))
})

test_that("posterior kernel computations work correctly", {
  kers <- posterior_predict_kernels(fit, X_pred)
  expect_equal(names(kers), c("data_vs_data", "pred_vs_data", "pred_vs_pred"))
  expect_equal(dim(kers$data_vs_data), c(5, 2, 24, 24))
  expect_equal(dim(kers$pred_vs_data), c(5, 2, 44, 24))
  expect_equal(dim(kers$pred_vs_pred), c(5, 2, 44, 44))
})
