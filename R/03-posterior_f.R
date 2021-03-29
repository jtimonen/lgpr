#' Evaluate or extract function posterior distributions
#'
#' @description Evaluates (or extracts) the posterior distribution of the
#' total signal \code{f} and its additive components (or draws from these
#' distributions). All these are computed for
#' each parameter draw (defined by \code{draws}), or other parameter set
#' (obtained by a reduction defined by \code{reduce}).
#'
#' @param fit An object of class \linkS4class{lgpfit}.
#' @param x A data frame of points where the predictions are computed.
#' The function \code{\link{new_x}} can help in creating it.
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
#' @return An object of class \linkS4class{FunctionPosterior} or
#' \linkS4class{FunctionDraws}.
#' @param verbose Should more some informational messages be printed?
#' @param debug_km Should this only return the required kernel matrices.
#' Can be used for debugging or testing \code{\link{kernels_fpost}}.
#' @param debug_dims Should this print dimensions of some variables.
#' Can be used for debugging \code{\link{kernels_fpost}}.
#' @param STREAM an external pointer
posterior_f <- function(fit,
                        x = NULL,
                        c_hat_pred = NULL,
                        reduce = function(x) base::mean(x),
                        draws = NULL,
                        verbose = TRUE,
                        debug_km = FALSE,
                        debug_dims = FALSE,
                        STREAM = get_stream()) {

  # Settings
  if (is.null(x)) x <- get_data(fit)
  if (!is.null(draws)) reduce <- NULL
  f_sampled <- is_f_sampled(fit)

  # Compute all required kernel matrices
  km <- kernels_fpost(fit, x, reduce, draws, verbose, debug_dims, STREAM)
  if (debug_km) {
    return(km)
  }

  # Compute the function posteriors
  if (f_sampled) {
    chp <- c_hat_pred
    out <- fp_latent(km, fit, x, chp, reduce, draws, verbose, STREAM)
  } else {
    out <- fp_marginal(km, fit, x, reduce, draws, verbose, STREAM)
  }
  return(out)
}

# Analytic function posteriors
fp_marginal <- function(km, fit, x, reduce, draws, verbose, STREAM) {

  # Helper function 1
  format_pred_comp <- function(m, s, comp_names) {
    S <- dim(m)[2] # number of parameter sets
    D <- dim(m)[3] # number of components
    P <- dim(m)[4] # number of prediction points
    paramset <- as.factor(rep(1:S, D * P))
    component <- as.factor(rep(rep(comp_names, each = S), P))
    eval_point <- as.factor(rep(1:P, each = D * S))
    m <- as.vector(m)
    s <- as.vector(s)
    check_lengths(m, s)
    check_lengths(m, paramset)
    check_lengths(m, component)
    check_lengths(m, eval_point)
    df <- data.frame(paramset, component, eval_point, m, s)
    colnames(df) <- c("paramset", "component", "eval_point", "mean", "std")
    return(df)
  }

  # Helper function 2
  format_pred_total <- function(m, s, sigma) {
    S <- dim(m)[2] # number of param sets
    P <- dim(m)[3] # number of prediction points
    paramset <- as.factor(rep(1:S, P))
    eval_point <- as.factor(rep(1:P, each = S))
    sigma <- rep(sigma, P)
    m <- as.vector(m)
    s <- as.vector(s)
    check_lengths(m, s)
    check_lengths(m, paramset)
    check_lengths(m, eval_point)
    check_lengths(m, sigma)
    df <- data.frame(paramset, eval_point, m, s, sigma)
    colnames(df) <- c("paramset", "eval_point", "mean", "std", "sigma")
    return(df)
  }

  # TODO: pred eqs
  fp <- km
  return(fp)

  # Format components
  # comp_names <- component_names(fit@model)
  # mc <- dollar(ext, "f_comp_mean")
  # sc <- dollar(ext, "f_comp_std")
  # df_comp <- format_pred_comp(mc, sc, comp_names)

  # Format total
  # m <- dollar(ext, "f_mean")
  # s <- dollar(ext, "f_std")
  ## sigma <- as.vector(dollar(stan_data, "d_sigma"))
  # df_total <- format_pred_total(m, s, sigma)

  # Return
  # new("FunctionPosterior",
  #  components = df_comp,
  #  total = df_total,
  #  x = x,
  #  model = fit@model
  # )
}

# Function posterior draws
fp_latent <- function(km, fit, x, c_hat_pred, reduce, draws, verbose, STREAM) {
  fpred <- NULL
  stop("FP_LATENT: NOT IMPLEMENTED!") # TODO

  # Create the prediction
  # h <- pred.latent_h(fit, f_pred, c_hat_pred, verbose)
  #
  ## Return
  # new("Prediction",
  #    f_comp = arr3_to_list(dollar(kr, "f_comp")),
  #    f = f,
  #    h = h
  # )
  return(f_pred)
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
