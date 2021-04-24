# A model with two age components
library(lgpr)
library(bench)

# Create data using two age components
setup_data <- function(N = 10, H = 10, sigma = 0.2) {
  ttt <- seq(12, 72, by = H)
  K <- length(ttt)
  age <- rep(ttt, N) + rep(H / N * 1:N, each = K)
  id <- as.factor(rep(1:N, each = K))
  dat <- data.frame(id, age)
  f1 <- exp(-sin(0.2 * age))
  f2 <- -0.001 * (-60 + age)**2
  f1 <- f1 / sd(f1)
  f2 <- f2 / sd(f2)
  dat$age2 <- age
  dat$y <- f1 + f2 + sigma * rnorm(n = length(age))
  return(dat)
}

# Create model where there are different priors for the two lengthscales
setup_model <- function(dat) {
  d1 <- log_normal(1.0, 0.3)
  d2 <- normal(0, 0.3)
  prior <- list(ell = list(d1, d2))
  model <- create_model(y ~ age + age2, dat, prior = prior)
  return(model)
}

# Run the experiment
run_test <- function(...) {
  dat <- setup_data()
  mod <- setup_model(dat)
  fit <- sample_model(mod, ...)
  return(fit)
}

# Run post-fitting tasks
run_post_tasks <- function(fit) {

  # Compute and plot predictions at data points
  p <- pred(fit, reduce = NULL)
  plt1 <- plot_components(fit, p, group_by = NA)
  plt2 <- plot_pred(fit, p, group_by = NA)

  # Compute and plot out-of-sample predictions
  x_pred <- new_x(fit, x_values = seq(0, 200, by = 1), group_by = NA)
  x_pred$age2 <- x_pred$age
  bm <- bench::mark(
    {
      pp <- pred(fit, x_pred, reduce = NULL)
    },
    iterations = 1
  )
  plt3 <- plot_components(fit, pp, group_by = NA)
  plt4 <- plot_pred(fit, pp, group_by = NA)

  # Return results
  results <- list(
    n_pred = nrow(x_pred),
    pred_time = as.character(bm$total_time),
    pred_mem = as.character(bm$mem_alloc),
    pred_bm = bm,
    plots = list(plt1, plt2, plt3, plt4)
  )
  return(results)
}
