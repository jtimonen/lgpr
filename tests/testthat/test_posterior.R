library(lgpr)

# -------------------------------------------------------------------------

context("Function posteriors and posterior predictive distributions")

my_prior <- list(
  alpha = normal(1, 0.1),
  ell = normal(1, 0.1),
  sigma = normal(1, 0.1)
)

dat <- testdata_001
model <- create_model(y ~ gp(age), dat, prior = my_prior)
fit <- sample_model(
  model = model,
  iter = 1200,
  chains = 1,
  refresh = 0,
  seed = 123
)

t <- seq(0, 50, by = 1)
X_pred <- new_x(dat, t)
a <- posterior_predict(fit, X_pred)

test_that("kernel parameter samples can be extracted", {
  a <- get_draws_kernel_pars(fit)
  expect_equal(names(a), c("alpha", "ell", "wrp", "beta", "teff"))
  expect_equal(nrow(a$alpha), nrow(a$teff))
})
