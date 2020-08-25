library(lgpr)
library(rstan)

# Create test data
N <- 6
a <- seq(10, 60, by = 10)
L <- length(a)
age <- rep(a, N)
id <- as.factor(rep(1:N, each = L))
sex <- as.factor(c(rep("Male", 3 * L), rep("Female", 3 * L)))
dis_age <- 32 - age
dis_age[id %in% c(1, 4, 5)] <- NA
y <- sin(10 * age) + 0.5 * as.numeric(sex)
dat <- data.frame(age, dis_age, id, sex, y)

# -------------------------------------------------------------------------

context("Main function lgp_fit")

test_that("stan model can be sampled with minimal input", {
  m <- lgp_model(y ~ gp(age) + categ(sex), dat)
  suppressWarnings({
    fit <- lgp_sampling(m, chains = 1, iter = 100, refresh = 0)
    lp <- rstan::extract(fit)$lp__
    expect_equal(length(lp), 50)
  })
})

test_that("stan model can be sampled with various components", {
  m <- lgp_model(y ~ heter(id) * gp(age) + gp_warp(dis_age), dat)
  stan_data <- m@stan_input

  suppressWarnings({
    fit <- lgp_sampling(m, chains = 1, iter = 10, refresh = 0)
    lp <- rstan::extract(fit)$lp__
    expect_equal(length(lp), 5)
  })
})
