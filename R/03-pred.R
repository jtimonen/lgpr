#' Compute posterior predictive distribution or get draws from it
#'
#' @export
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param x A data frame of points where the predictions are computed.
#' The function \code{\link{new_x}} provides an easy way to create it.
#' @param c_hat_pred This is only used if the latent signal \code{f} was
#' sampled. This input contains the values added to the sum \code{f} before
#' passing through inverse link function. Must be a vector with length equal to
#' the number of prediction points. If original \code{c_hat} was constant,
#' then \code{c_hat_pred} can be ignored, in which case this will by default
#' use the same constant.
#' @param reduce Reduction for parameters draws. Can be a function that
#' is applied to reduce all parameter draws into one parameter set, or
#' NULL (no reduction). Has no effect if \code{draws} is specified.
#' @param draws Indices of parameter draws to use, or \code{NULL} to use all
#' draws.
#' @param STREAM External pointer. By default obtained with
#' \code{\link[rstan]{get_stream}}.
#' @return An object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}.
#' @param ... Other optional arguments passed to \code{\link{posterior_f}}.
#' @param verbose Should more information and a progress bar be printed?
pred <- function(fit, x, c_hat_pred = NULL,
                 reduce = function(x) base::mean(x), draws = NULL,
                 verbose = TRUE, STREAM = get_stream(), ...) {
  check_type(x, "data.frame")
  f_sampled <- is_f_sampled(fit)
  if (class(fit) == "lgpfit") {
    fp <- posterior_f(fit, x, reduce, draws, verbose, STREAM, ...)
    return(fp)
  }
  return("not implemented!")
}


# Transform distribution of f to distribution of y
pred_marginal.f_to_y <- function(f_mean, f_std, sigma, y_norm_inv) {
  y_mean <- f_mean
  y_var <- add_to_columns(f_std^2, sigma^2)
  y_std <- sqrt(y_var)
  y_upper <- y_mean + y_std

  # Scale y_pred to original scale
  y_mean <- call_fun(y_norm_inv, y_mean)
  y_upper <- call_fun(y_norm_inv, y_upper)
  y_std <- y_upper - y_mean

  # Return
  list(mean = y_mean, std = y_std)
}

# Map the sum f from pred.latent_compute to h
pred.latent_h <- function(fit, f, c_hat_pred, verbose) {

  # helper function
  is_constant <- function(x) {
    s <- sum(x == x[1])
    s == length(x)
  }

  num_draws <- dim(f)[1]
  num_pred <- dim(f)[2]
  if (is.null(c_hat_pred)) {
    c_hat_data <- get_chat(fit)
    if (is_constant(c_hat_data)) {
      msg <- paste0(
        "c_hat_pred not given,",
        " using constant c_hat_pred = ", c_hat_data[1], "\n"
      )
      if (verbose) cat(msg)
      c_hat_pred <- rep(c_hat_data[1], num_pred)
    } else {
      msg <- paste0(
        "c_hat (at data points) is not constant! ",
        "you must give c_hat_pred (at prediction points) as input!"
      )
      stop(msg)
    }
  }

  f <- f + repvec(c_hat_pred, num_draws)
  h <- link_inv(f, get_obs_model(fit))
  return(h)
}
