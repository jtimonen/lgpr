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
DRAW_INDS <- round(seq(1, NUM_ITER * NUM_CHAINS / 2, length.out = 100))
verbose <- FALSE

# Run test suite
suite_path <- file.path("tests", "suite")
source(file.path(suite_path, "common.R"))

examples_path <- file.path(suite_path, "examples")
files <- dir(examples_path)
INFO <- c()

for (f in files) {

  # Setup
  fp <- file.path(examples_path, f)
  cat("Running:", f, "\n")
  source(fp)
  start_time <- Sys.time()

  # Run model fitting
  res_fit <- run_example(
    verbose,
    iter = NUM_ITER, chains = NUM_CHAINS, cores = NUM_CORES,
    refresh = REFRESH, seed = STAN_SEED
  )
  fit <- res_fit$fit
  t_fit <- res_fit$time

  # Run post-fitting tasks
  t_pred <- run_pred(fit, verbose)
  t_post <- run_post_tasks(fit, DRAW_INDS, verbose)

  # Store info
  t_total <- as.numeric(Sys.time() - start_time)
  info <- get_info(fit, f, t_fit, t_pred, t_post, t_total)
  INFO <- rbind(INFO, info)
}
