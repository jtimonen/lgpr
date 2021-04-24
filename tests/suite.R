# A larger test suite
#  - tests correctness of inference
#  - measures runtime
library(profvis)
library(lgpr)
library(rstan)

NUM_ITER <- 2000
NUM_CHAINS <- 4
NUM_CORES <- 4
REFRESH <- 0
STAN_SEED <- 123
DRAW_INDS <- round(seq(1, NUM_ITER*NUM_CHAINS/2, length.out = 100))

# Get experiment information from lgpfit object
get_info <- function(fit, name, n_pred, pred_time1, pred_time2, total_time) {
  n_obs <- lgpr:::get_num_obs(fit)
  n_comps <- lgpr:::get_num_comps(fit)
  f_sampled <- is_f_sampled(fit)
  n_div <- rstan::get_num_divergent(fit@stan_fit)
  chain_times <- apply(rstan::get_elapsed_time(fit@stan_fit), 1, mean)
  fit_time_mean <- mean(chain_times)
  fit_time_sd <- stats::sd(chain_times)
  df <- data.frame(
    name, f_sampled, n_div, n_obs, n_comps, fit_time_mean, fit_time_sd,
    n_pred, pred_time1, pred_time2, total_time
  )
  return(df)
}

# Run test suite
suite_path <- file.path("tests/suite")
files <- dir(suite_path)
INFO <- c()

for (f in files) {

  # Setup
  fp <- file.path(suite_path, f)
  cat("FILE:", fp, "\n")
  source(fp)
  start_time <- Sys.time()

  # Run model fitting
  fit <- run_test(
    iter = NUM_ITER, chains = NUM_CHAINS, cores = NUM_CORES,
    refresh = REFRESH, seed = STAN_SEED
  )

  # Run post-fitting tasks
  res <- run_post_tasks(fit, DRAW_INDS)

  # Store info
  total_time <- as.numeric(Sys.time() - start_time)
  info <- get_info(
    fit, f, res$n_pred, res$time1, res$time2,
    total_time
  )
  INFO <- rbind(INFO, info)
}
