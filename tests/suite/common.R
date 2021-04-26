# Object size
r_size_kb <- function(object) {
  format(object.size(object), units = "Kb")
}

# File size
file_size_kb <- function(file) {
  fi <- file.info(file)
  fs_mb <- fi$size / (100)
  paste0(round(fs_mb, 1), " Kb")
}

# Get experiment information from lgpfit object
get_info <- function(fit, name, t_fit, t_pred, t_post, t_total, size_disk) {
  n_obs <- lgpr:::get_num_obs(fit)
  n_comps <- lgpr:::get_num_comps(fit)
  obs_model <- lgpr:::get_obs_model(fit)
  f_sampled <- is_f_sampled(fit)
  n_div <- rstan::get_num_divergent(fit@stan_fit)
  mf <- rstan::monitor(fit@stan_fit, print = FALSE)
  min_TESS <- min(mf$Tail_ESS)
  min_BESS <- min(mf$Bulk_ESS)
  max_Rhat <- max(mf$Rhat)
  chain_times <- apply(rstan::get_elapsed_time(fit@stan_fit), 1, mean)
  t_ch_mean <- mean(chain_times)
  t_ch_sd <- stats::sd(chain_times)
  size_fit <- r_size_kb(fit)
  size_small <- r_size_kb(clear_postproc(fit))

  # Return data frame
  data.frame(
    name, obs_model, f_sampled, n_obs, n_comps,
    t_ch_mean, t_ch_sd, t_fit, t_pred, t_post, t_total,
    size_fit, size_small, size_disk,
    n_div, max_Rhat, min_TESS, min_BESS
  )
}

# Run an example
run_example <- function(verbose, ...) {
  mod <- setup_model()
  start_time <- Sys.time()
  fit <- sample_model(mod, quiet = !verbose, ...)
  elapsed_time <- as.double(Sys.time() - start_time, units = "secs")
  list(fit = fit, time = elapsed_time)
}

# Time pred() at data points with reduce = NULL
run_pred <- function(fit, verbose) {
  x_pred <- get_data(fit)
  start_time <- Sys.time()
  p <- pred(fit, x = x_pred, reduce = NULL, verbose = verbose)
  elapsed_time <- as.double(Sys.time() - start_time, units = "secs")
  return(elapsed_time)
}

# Format result data frame
round_results <- function(info, t_digits, rhat_digits) {
  t_cols <- c("t_ch_mean", "t_ch_sd", "t_fit", "t_pred", "t_post", "t_total")
  info[t_cols] <- round(INFO[t_cols], t_digits)
  info["max_Rhat"] <- round(info["max_Rhat"], rhat_digits)
  return(info)
}
