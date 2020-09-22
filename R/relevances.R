#' Assess component relevances
#'
#' @description
#' \itemize{
#'   \item \code{relevances} returns a named vector with length equal to
#'   number of components plus one
#'   \item \code{relevances.default} is the default method
#'   \item \code{relevances.default.all} is a helper function
#' }
#' @param fit an object of class \code{lgpfit}
#' @param ... currently has no effect
#' @name relevances
NULL

#' @export
#' @rdname relevances
relevances <- function(fit, ...) {
  check_type(fit, "lgpfit")
  relevances.default(fit, ...)
}

#' @rdname relevances
relevances.default <- function(fit, ...) {
  rels <- relevances.default.all(fit)
  colMeans(rels)
}

#' @rdname relevances
relevances.default.all <- function(fit) {
  pred <- get_pred(fit)
  p_noise <- relevances.default.noise(fit, pred)
  p_comps <- relevances.default.comp(fit, pred)
  num_comps <- ncol(p_comps)
  mult <- t(repvec(1 - p_noise, num_comps))
  p_comps <- mult * p_comps
  df <- cbind(p_comps, p_noise)
  colnames(df) <- c(names(p_comps), "noise")
  return(df)
}

#' Helper functions for relevances.default
#'
#' @description
#' \itemize{
#'   \item \code{relevances.default.noise} computes the noise proportion
#'   in each draws, and returns a vector with length \code{num_draws}
#'   \item \code{relevances.default.comp} divides the signal proportion for
#'   each component in each draw, and returns a data frame with
#'   \code{num_draws} rows and \code{num_comps} colums
#' }
#' @param pred an object of class \linkS4class{Prediction} or
#' \linkS4class{GaussianPrediction}
#' @inheritParams relevances
#' @name relevances_default
#' @return an array or vector
NULL

#' @rdname relevances_default
relevances.default.noise <- function(fit, pred) {
  typ <- class(pred)
  h <- if (typ == "Prediction") pred@h else pred@y_mean
  num_obs <- ncol(h)
  num_draws <- nrow(h)
  y <- get_y(fit, original = TRUE)
  resid <- h - repvec(y, num_draws)
  signal_var <- row_vars(h)
  error_var <- rowSums(resid^2) / (num_obs - 1)
  p_signal <- signal_var / (signal_var + error_var)
  return(1 - p_signal)
}

#' @rdname relevances_default
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
