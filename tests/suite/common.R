# Get experiment information from lgpfit object
get_info <- function(fit, name, t_fit, t_pred, t_post, t_total) {
  n_obs <- lgpr:::get_num_obs(fit)
  n_comps <- lgpr:::get_num_comps(fit)
  f_sampled <- is_f_sampled(fit)
  n_div <- rstan::get_num_divergent(fit@stan_fit)
  chain_times <- apply(rstan::get_elapsed_time(fit@stan_fit), 1, mean)
  t_ch_mean <- mean(chain_times)
  t_ch_sd <- stats::sd(chain_times)
  df <- data.frame(
    name, f_sampled, n_div, n_obs, n_comps, t_ch_mean, t_ch_sd,
    t_fit, t_pred, t_post, t_total
  )
  return(df)
}

# Run an example
run_example <- function(verbose, ...) {
  dat <- setup_data()
  mod <- setup_model(dat)
  start_time <- Sys.time()
  fit <- sample_model(mod, quiet = !verbose, ...)
  elapsed_time <- as.numeric(Sys.time() - start_time)
  list(fit = fit, time = elapsed_time)
}

# Time pred() at data points with reduce = NULL
run_pred <- function(fit, verbose) {
  x_pred <- get_data(fit)
  start_time <- Sys.time()
  p <- pred(fit, x = x_pred, reduce = NULL, verbose = verbose)
  elapsed_time <- as.numeric(Sys.time() - start_time)
  return(elapsed_time)
}
