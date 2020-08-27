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

context("Sampling or optimizing a model")

test_that("lgpmodel can be sampled with minimal input", {
  m <- create_model(y ~ gp(age) + categ(sex), dat)
  suppressWarnings({
    fit <- sample_model(m, chains = 1, iter = 100, refresh = 0)
    lp <- rstan::extract(fit)$lp__
    expect_equal(length(lp), 50)
    expect_s4_class(fit, "stanfit")
  })
})

test_that("lgpmodel with various components can be sampled", {
  m <- create_model(y ~ heter(id) * gp(age) + gp_warp(dis_age), dat)

  suppressWarnings({
    fit <- sample_model(m, chains = 1, iter = 10, refresh = 0)
    lp <- rstan::extract(fit)$lp__
    expect_equal(length(lp), 5)
  })
})


test_that("lgpmodel can be optimized", {
  m <- create_model(y ~ heter(id) * gp(age), dat)

  suppressWarnings({
    fit <- optimize_model(m, iter = 10)
    # names(fit) = c("par", "value", "return_code", "theta_tilde")
    expect_equal(class(fit), "list")
  })
})
