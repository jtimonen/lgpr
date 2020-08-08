library(lgpr)
library(rstan)
sm <- lgpr::get_stan_model()

context("Main stan model")

test_that("stan model can be sampled with minimal input", {
  stan_data <- lgpr:::minimal_stan_data()
  suppressWarnings({
    fit <- rstan::sampling(sm,
      data = stan_data, chains = 1,
      iter = 100, refresh = 0
    )
    lp <- extract(fit)$lp__
    expect_equal(length(lp), 50)
  })
})
