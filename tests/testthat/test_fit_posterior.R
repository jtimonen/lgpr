library(lgpr)

# -------------------------------------------------------------------------

context("Posterior sampling and optimization")

my_prior <- list(
  alpha = normal(1, 0.1),
  ell = normal(1, 0.1),
  sigma = normal(1, 0.1)
)

model <- create_model(y ~ gp(age), testdata_001, prior = my_prior)

test_that("sample_model can do posterior sampling", {
  fit <- sample_model(
    model = model,
    iter = 1200,
    chains = 1,
    refresh = 0,
    seed = 123
  )

  expect_s4_class(fit, "lgpfit")
  p1 <- plot_draws(fit)
  p2 <- plot_draws(fit, type = "trace")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("optimize_model can optimize MAP parameters", {
  fit <- optimize_model(model, iter = 10, seed = 123)
  expect_equal(class(fit), "list")
})
