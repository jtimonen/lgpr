# A model with two age components
library(lgpr)

# Create data using two age components
setup_data <- function(N = 10, H = 10, sigma = 0.2) {
  ttt <- seq(12, 72, by = H)
  K <- length(ttt)
  age <- rep(ttt, N) + rep(H / N * 1:N, each = K)
  id <- as.factor(rep(1:N, each = K))
  dat <- data.frame(id, age)
  f1 <- -0.001 * (-60 + age)**2 # slowly changing
  f2 <- exp(-sin(0.3 * age)) # quickly changing
  f1 <- f1 / sd(f1)
  f2 <- f2 / sd(f2)
  dat$age_fast <- age
  dat$y <- f1 + f2 + sigma * rnorm(n = length(age))
  return(dat)
}

# Create model where there are different priors for the two lengthscales
setup_model <- function(dat) {
  d1 <- log_normal(1.0, 0.3)
  d2 <- normal(0, 0.3)
  prior <- list(ell = list(d1, d2))
  model <- create_model(y ~ age + age_fast, dat, prior = prior)
  return(model)
}

# Run post-fitting tasks
run_post_tasks <- function(fit, draw_inds, verbose) {
  start_time <- Sys.time()

  # Compute and plot predictions at data points
  p <- pred(fit, reduce = NULL, draws = draw_inds, verbose = verbose)
  plt1 <- plot_components(fit, p, group_by = NA, draw = FALSE)
  plt2 <- plot_pred(fit, p, group_by = NA)

  # Compute and plot out-of-sample predictions
  x_pred <- new_x(fit, x_values = seq(0, 200, by = 1), group_by = NA)
  x_pred$age_fast <- x_pred$age

  pp <- pred(fit, x_pred, reduce = NULL, draws = draw_inds, verbose = verbose)

  plt3 <- plot_components(fit, pp, group_by = NA, draw = FALSE)
  plt4 <- plot_pred(fit, pp, group_by = NA)

  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}

# Run relevance and selection tasks
test_relevances_and_selection <- function(fit) {
  r <- relevances(fit)
  s <- select(fit)
  print(r)
  print(s)
}
