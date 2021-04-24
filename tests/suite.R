# A larger test suite
#  - tests correctness of inference
#  - measures runtime
library(bench)
library(lgpr)
library(rstan)

NUM_ITER <- 2000
NUM_CHAINS <- 4
NUM_CORES <- 4
REFRESH <- 0

# Get experiment information from lgpfit object
get_info <- function(fit, iter, chains, cores, name, total_time, mem_alloc) {
  obs <- lgpr:::get_num_obs(fit)
  comps <- lgpr:::get_num_comps(fit)
  f_sampled <- is_f_sampled(fit)
  div <- rstan::get_num_divergent(fit@stan_fit)
  times <- apply(rstan::get_elapsed_time(fit@stan_fit), 1, mean)
  t_mean <- mean(times)
  t_sd <- stats::sd(times)
  df <- data.frame(
    name, iter, chains, cores, div, f_sampled,
    obs, comps, t_mean, t_sd, total_time, mem_alloc
  )
  return(df)
}

# Run test suite
suite_path <- file.path("tests/suite")
files <- dir(suite_path)
for (f in files) {
  fp <- file.path(suite_path, f)
  cat("FILE:", fp, "\n")
  source(fp)
  fit <- run_test(
    iter = NUM_ITER, chains = NUM_CHAINS, cores = NUM_CORES,
    refresh = REFRESH
  )
  bm <- bench::mark(
    {
      r <- relevances(fit)
    },
    iterations = 1
  )
  mem <- bm$mem_alloc
  tt <- bm$total_time
  info <- get_info(fit, NUM_ITER, NUM_CHAINS, NUM_CORES, f, tt, mem)
}
