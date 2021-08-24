#' Posterior predictions and function posteriors
#'
#' @description
#' \itemize{
#'   \item If \code{use_f_draws} is \code{FALSE}, this
#'   computes the conditional posterior
#'   distributions of each model component, their sum, and the conditional
#'   predictive distribution. All these are computed given
#'   each (hyper)parameter draw (defined by \code{draws}), or other parameter
#'   set (obtained by a reduction defined by \code{reduce}). Results are stored
#'   in a \linkS4class{GaussianPrediction} object which is then returned.
#'
#'   \item If \code{use_f_draws} is \code{TRUE}, this will
#'   extract the function component samples, compute their sum, and a
#'   version of the sum \code{f} that is transformed through the inverse
#'   link function. If \code{x} is not \code{NULL}, the function draws are
#'   extrapolated to the points specified by \code{x} using kernel regression.
#'   Results are stored in a \linkS4class{Prediction}
#'   object which is then returned.
#' }
#'
#' @export
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param x A data frame of points where function posterior distributions
#' and predictions should be computed or sampled.
#' The function \code{\link{new_x}} provides an easy way to create it.
#' If this is \code{NULL}, the data points are used.
#' @param c_hat_pred This is only used if the latent signal \code{f} was
#' sampled. This input contains the values added to the sum \code{f} before
#' passing through inverse link function. Must be a vector with length equal to
#' the number of prediction points. If original \code{c_hat} was constant,
#' then \code{c_hat_pred} can be ignored, in which case this will by default
#' use the same constant.
#' @param reduce Reduction for parameters draws. Can be a function that
#' is applied to reduce all parameter draws into one parameter set, or
#' \code{NULL} (no reduction). Has no effect if \code{draws} is specified.
#' @param draws Indices of parameter draws to use, or \code{NULL} to use all
#' draws.
#' @param STREAM External pointer. By default obtained with
#' \code{rstan::get_stream()}.
#' @param verbose Should more information and a possible progress bar be
#' printed?
#' @param use_f_draws Should function posterior draws be used. These only
#' exist if they were sampled in the first place (i.e. f was not marginalized).
#' @param force This is by default \code{FALSE} to prevent unintended
#' large computations that might crash R or take forever. Set it to \code{TRUE}
#' try computing anyway.
#' @return An object of class \linkS4class{GaussianPrediction} or
#' \linkS4class{Prediction}.
#' @family main functions
pred <- function(fit,
                 x = NULL,
                 reduce = function(x) base::mean(x),
                 draws = NULL,
                 use_f_draws = is_f_sampled(fit),
                 verbose = TRUE,
                 STREAM = get_stream(),
                 c_hat_pred = NULL,
                 force = FALSE) {
  vrb <- verbose
  if (!is.null(draws)) reduce <- NULL
  is_marg <- is_f_marginalized(fit)
  if (is_marg && use_f_draws) {
    stop("f was not sampled, so cannot set use_f_draws=TRUE!")
  }

  # If x=NULL, return existing pred or use x = data
  check_type(fit, "lgpfit")
  if (is.null(x)) {
    if (contains_postproc(fit)) {
      log_info("Using existing postprocessing information.", verbose)
      pr <- fit@postproc_results
      return(dollar(pr, "pred"))
    } else {
      log_info("Using x = data since x was NULL.", verbose)
      x <- get_data(fit)
      x_is_data <- TRUE
    }
  } else {
    prevent_too_large_mats(fit, x, reduce, draws, vrb, force)
    x_is_data <- FALSE
  }
  check_not_null(x)

  # Call subroutine
  if (!use_f_draws) {
    if (get_obs_model(fit) != "gaussian") {
      stop("Cannot set use_f_draws=FALSE unless likelihood is gaussian!")
    }
    fp <- posterior_f_gaussian(fit, x, x_is_data, reduce, draws, vrb, STREAM)
    out <- posterior_pred_gaussian(fit, fp, verbose)
  } else {
    fp <- posterior_f_latent(fit, x, x_is_data, reduce, draws, vrb, STREAM)
    out <- posterior_pred_latent(fit, fp, c_hat_pred, verbose)
  }
  return(out)
}


# Gaussian function posteriors to posterior pred
posterior_pred_gaussian <- function(fit, fp, verbose) {
  f_mean <- dollar(fp, "f_mean")
  f_std <- dollar(fp, "f_std")
  sigma2 <- dollar(fp, "sigma2")
  y_scl <- fit@model@y_scaling
  y_pred <- map_f_to_y(f_mean, f_std, sigma2, y_scl)
  new("GaussianPrediction",
    f_comp_mean = dollar(fp, "f_comp_mean"),
    f_comp_std = dollar(fp, "f_comp_std"),
    f_mean = f_mean,
    f_std = f_std,
    y_mean = dollar(y_pred, "mean"),
    y_std = dollar(y_pred, "std"),
    x = dollar(fp, "x")
  )
}

# pred when sample_f = TRUE
posterior_pred_latent <- function(fit, fp, c_hat_pred, verbose) {
  f <- dollar(fp, "sum")
  model <- fit@model
  c_hat_pred <- set_c_hat_pred(model, f, c_hat_pred, verbose)
  h <- map_f_to_h(fit, f, c_hat_pred, reduce = NULL)

  # Return
  new("Prediction",
    f_comp = dollar(fp, "comp"),
    f = f,
    h = h,
    x = dollar(fp, "x"),
    extrapolated = dollar(fp, "extrapolated")
  )
}

# Set c_hat_pred
set_c_hat_pred <- function(model, f, c_hat_pred, verbose) {

  # helper function
  is_constant <- function(x) {
    s <- sum(x == x[1])
    s == length(x)
  }

  num_pred <- dim(f)[2]
  if (is.null(c_hat_pred)) {
    c_hat_data <- get_chat(model)
    if (is_constant(c_hat_data)) {
      msg <- paste0(
        "c_hat_pred not given,",
        " using constant c_hat_pred = ", c_hat_data[1], "\n"
      )
      log_info(msg, verbose)
      c_hat_pred <- rep(c_hat_data[1], num_pred)
    } else {
      msg <- paste0(
        "c_hat (at data points) is not constant! ",
        "you must give c_hat_pred (at prediction points) as input!"
      )
      stop(msg)
    }
  }
  return(c_hat_pred)
}
