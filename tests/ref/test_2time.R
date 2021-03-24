library(lgpr)

context("Full analysis using two kernels for age and group_by = NA")


a <- simulate_data(
  N = 8,
  t_data = seq(12, 96, by = 12),
  covariates = c(1, 1, 1),
  t_jitter = 4
)
dat <- a@data

age <- dat$age
f1 <- exp(-sin(0.2 * age))
f2 <- -0.001 * (-60 + age)**2

dat$age2 <- age
dat$y <- f1 + f2 + 0.5 * rnorm(n = length(age))

# Create model where there are different priors for the two engthscales
d1 <- normal(1, 0.3)
d2 <- normal(0, 0.2)
prior <- list(ell = list(d1, d2))
model <- create_model(y ~ age + age2, dat, prior = prior, sample_f = TRUE)

test_that("different priors can be set for parameters of same type", {
  hp_ell <- model@stan_input$hyper_ell
  expect_equal(hp_ell[1, 1], 1.0)
  expect_equal(hp_ell[1, 2], 0.3)
  expect_equal(hp_ell[2, 1], 0.0)
  expect_equal(hp_ell[2, 2], 0.2)
})

suppressWarnings({
  fit <- sample_model(model, iter = 10, chains = 1, refresh = 0)
})

test_that("prediction and plotting is possible with group_by = NA", {
  p <- get_pred(fit, reduce = mean)
  expect_s4_class(p, "Prediction")

  plt1 <- plot_components(fit, group_by = NA, reduce = NULL)
  plt2 <- plot_pred(fit, group_by = NA, reduce = NULL)
  expect_equal(length(plt1), 3)
  expect_s3_class(plt2, "ggplot")

  x_pred <- new_x(dat, x_values = seq(0, 200, by = 2), group_by = NA)
  x_pred$age2 <- x_pred$age
  pp <- pred(fit, x_pred, reduce = NULL, draws = c(1, 2, 5), verbose = FALSE)

  plt3 <- plot_components(fit, x_pred, pp, group_by = NA)
  plt4 <- plot_pred(fit, x_pred, pp, group_by = NA)
  expect_equal(length(plt3), 3)
  expect_s3_class(plt4, "ggplot")
})
