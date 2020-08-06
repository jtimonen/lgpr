library(lgpr)
library(rstan)
sm <- lgpr::get_stan_model()

test_that("stan model can be sampled with minimal input", {
  stan_data <- lgpr:::minimal_stan_data()
  print(stan_data)
  fit <- rstan::sampling(sm, data = stan_data, chains = 1, iter = 100)
  expect_equal(1, 1)
})

