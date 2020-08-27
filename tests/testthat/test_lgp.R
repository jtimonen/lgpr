library(lgpr)

# -------------------------------------------------------------------------

context("Main function lgp")

# Create test input
sim <- simulate_data(
  N = 4,
  t_data = seq(6, 36, by = 6),
  covariates = c(0, 1, 2, 3),
  lengthscales = rep(12, 5),
  relevances = rep(1, 6),
  t_jitter = 0.5
)

# Formula
formula <- y ~ zerosum(id) * gp(age) + gp_warp_vm(diseaseAge) +
  categ(z) + gp(age) + gp(x)

test_that("an lgpfit object is returned and can be plotted", {
  suppressWarnings({
    fit <- lgp(
      formula = formula,
      data = sim@data,
      iter = 100,
      chains = 2,
      refresh = 0, cores = 2
    )
    expect_s4_class(fit, "lgpfit")
    p1a <- plot_posterior(fit)
    p1b <- plot_posterior(fit, type = "trace")
    p1c <- plot_posterior(fit, type = "dens", regex_pars = "alpha")
    p2 <- plot_posterior_warp(fit)
    expect_s3_class(p1a, "ggplot")
    expect_s3_class(p1b, "ggplot")
    expect_s3_class(p1c, "ggplot")
    expect_s3_class(p2, "ggplot")
  })
})
