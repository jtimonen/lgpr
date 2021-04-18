#' Assess component relevances
#'
#' @export
#' @param fit an object of class \code{lgpfit}
#' @param reduce a function to apply to reduce the relevances given each
#' parameter draw into one value
#' @param verbose Can this print any messages?
#' @return a named vector with length equal to \code{num_comps + 1}
#' @param ... currently has no effect
relevances <- function(fit,
                       reduce = function(x) base::mean(x),
                       verbose = TRUE,
                       ...) {
  check_type(fit, "lgpfit")
  relevances.default(fit, reduce, verbose, ...)
}

# Default method
relevances.default <- function(fit, reduce, verbose, ...) {
  df <- relevances.default.all(fit, verbose)
  if (is.null(reduce)) {
    return(df)
  } else {
    check_type(reduce, "function")
    df <- apply(df, 2, reduce)
  }
  return(df)
}

# returns a data frame of size (num_draws) x (num_comps + 1)
relevances.default.all <- function(fit, verbose) {
  pred <- get_pred(fit, reduce = NULL, verbose = verbose)
  p_noise <- relevances.default.noise(fit, pred)
  p_comps <- relevances.default.comp(fit, pred)
  num_comps <- ncol(p_comps)
  mult <- t(repvec(1 - p_noise, num_comps))
  p_comps <- mult * p_comps
  df <- cbind(p_comps, p_noise)
  colnames(df) <- c(names(p_comps), "noise")
  return(df)
}

# Computes the noise proportion in each draw, and returns a vector with
# length num_draws
relevances.default.noise <- function(fit, pred) {
  typ <- class(pred)
  h <- if (typ == "Prediction") pred@h else pred@y_mean
  num_obs <- ncol(h)
  num_draws <- nrow(h)
  y <- get_y(fit, original = TRUE)
  y <- divide_by_num_trials(y, fit)
  resid <- h - repvec(y, num_draws)
  signal_var <- row_vars(h)
  error_var <- rowSums(resid^2) / (num_obs - 1)
  p_signal <- signal_var / (signal_var + error_var)
  return(1 - p_signal)
}

# Divides the signal proportion for each component in each draw, and returns
# a data frame with \code{num_draws} rows and \code{num_comps} colums
relevances.default.comp <- function(fit, pred) {
  typ <- class(pred)
  f_comp <- if (typ == "Prediction") pred@f_comp else pred@f_comp_mean
  v_list <- lapply(f_comp, row_vars)
  nam <- names(v_list)
  df <- data.frame(v_list)
  colnames(df) <- nam
  df <- normalize_rows(df)
  return(df)
}
