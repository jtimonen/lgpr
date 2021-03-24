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
#' NULL (no reduction). Has no effect if \code{draws} is specified.
#' @param draws Indices of parameter draws to use, or \code{NULL} to use all
#' draws.
#' @return An object of class \linkS4class{FunctionPosterior} or
#' \linkS4class{FunctionDraws}.
#' @param verbose Should more information be printed?
#' @param refresh How often to print progress? Has no effect if \code{verbose}
#' is \code{FALSE}.
posterior_f <- function(fit,
                        x = NULL,
                        c_hat_pred = NULL,
                        reduce = function(x) base::mean(x),
                        draws = NULL,
                        verbose = TRUE,
                        refresh = NULL) {
  if (is.null(x)) x <- get_data(fit)
  f_sampled <- is_f_sampled(fit)
  if (!verbose) refresh <- 0
  if (!is.null(draws)) reduce <- NULL
  if (f_sampled) {
    out <- fp_latent(fit, x, c_hat_pred, reduce, draws, refresh)
  } else {
    out <- fp_marginal(fit, x, reduce, draws, refresh)
  }
  return(out)
}

# Analytic function posteriors
fp_marginal <- function(fit, x, reduce, draws, refresh) {

  # Compute f_pred
  stan_data <- fp_input(fit, x, reduce, draws, refresh)
  stan_model <- dollar(stanmodels, "fp_marginal")
  stan_fit <- rstan::sampling(
    object = stan_model,
    data = stan_data,
    check_data = TRUE,
    algorithm = "Fixed_param",
    iter = 1,
    chains = 1,
    refresh = dollar(stan_data, "refresh")
  )

  # Helper function 1
  format_pred_comp <- function(m, s, comp_names) {
    S <- dim(m)[2] # number of param sets
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

  # Extract
  pars <- c("f_comp_mean", "f_mean", "f_comp_std", "f_std")
  ext <- rstan::extract(stan_fit, pars = pars)

  # Format components
  comp_names <- component_names(fit@model)
  mc <- dollar(ext, "f_comp_mean")
  sc <- dollar(ext, "f_comp_std")
  df_comp <- format_pred_comp(mc, sc, comp_names)

  # Format total
  m <- dollar(ext, "f_mean")
  s <- dollar(ext, "f_std")
  sigma <- as.vector(dollar(stan_data, "d_sigma"))
  df_total <- format_pred_total(m, s, sigma)

  # Return
  new("FunctionPosterior",
    components = df_comp,
    total = df_total,
    x = x,
    model = fit@model
  )
}

# Function posterior draws
fp_latent <- function(fit, x, c_hat_pred, reduce, draws, refresh) {

  # Compute f_pred which has dim = c(S, num_pred*num_comps)
  if (is.null(x)) {
    f_pred <- get_draws(fit, draws, reduce, pars = "f_latent")
  } else {
    stan_data <- fp_input(fit, x, reduce, draws, refresh)
    stan_model <- dollar(stanmodels, "fp_latent")
    stan_fit <- rstan::sampling(
      object = stan_model,
      data = stan_data,
      check_data = TRUE,
      algorithm = "Fixed_param",
      iter = 1,
      chains = 1,
      refresh = 0
    )
    f_pred <- rstan::extract(stan_fit, pars = "f_draws")$f_draws[1, 1, , ]
  }

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
