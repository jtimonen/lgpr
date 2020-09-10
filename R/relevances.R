#' Assess component relevances
#'
#' \itemize{
#'   \item \code{relevances} assesses covariate relevances
#'   \item \code{relevances.f_marginalized} assesses covariate relevances for a
#'   model where the latent function \code{f} was not sampled
#'   \item \code{relevances.f_sampled} assesses covariate relevances for a
#'   model where the latent function \code{f} was sampled
#' }
#' @param fit an object of class \code{lgpfit}
#' @name relevances
NULL

#' @rdname relevances
relevances <- function(fit) {
  check_type(fit, "lgpfit")
  flag <- is_f_sampled(fit)
  if (flag) {
    out <- relevances.f_sampled(fit)
  } else {
    out <- relevances.f_marginalized(fit)
  }
  return(out)
}

#' @rdname relevances
relevances.f_sampled <- function(fit) {
  get_component_info(fit)
}

#' @rdname relevances
relevances.f_marginalized <- function(fit) {
  y <- get_y(fit, original = FALSE)
  f <- get_f(fit)
  f_total_mean <- dollar(dollar(dollar(f, "f"), "total"), "mean")
  obs_model <- get_obs_model(fit)
  p_noise <- noise_proportion(f_total_mean, y, obs_model)
  p_noise
}


#' Determine noise proportion for each sample
#'
#' @param f_total an array of shape \code{num_samples} x \code{num_obs}
#' @param y an array of shape \code{1} x \code{num_obs}
#' @param likelihood used likelihood model name
#' @return a value between 0 and 1
noise_proportion <- function(f_total, y, likelihood) {
  num_draws <- nrow(f_total)
  num_obs <- ncol(f_total)
  g_total <- link_inv(f_total, likelihood)
  resid <- g_total - repvec(y, num_draws)

  signal_var <- apply(g_total, 1, stats::var)
  error_var <- colSums(resid^2) / (num_obs - 1)
  p_signal <- signal_var / (signal_var + error_var)
  return(1 - p_signal)
}
