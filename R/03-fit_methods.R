#' Graphical posterior or prior predictive checks
#'
#' @export
#' @inheritParams get_draws
#' @param data the original data frame
#' @param fun \code{bayesplot} function name
#' @param ... additional arguments passed to the default
#' \code{\link[bayesplot]{pp_check}} method in
#' \code{bayesplot}
#' @return a \code{ggplot} object
#' @seealso Introduction to graphical posterior predictive checks:
#' \href{here}{https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html}
#'
ppc <- function(fit, data, fun = default_ppc_fun(fit), ...) {
  check_type(fit, "lgpfit")
  check_type(data, "data.frame")
  check_type(fun, "function")
  y_name <- get_y_name(fit)
  y <- dollar(data, y_name)
  y_rep <- get_y_rng(fit, original_scale = TRUE)
  bayesplot::pp_check(y, y_rep, fun, ...)
}

#' Extract posterior or prior predictive distribution draws
#'
#' @export
#' @inheritParams get_draws
#' @param original_scale should the draws be scaled back to original y scale
#' (only has effect if likelihood is "gaussian", when data has been normalized
#' to zero mean and unit variance)
#' @return an array of shape \code{num_draws} x \code{num_obs}
#' @family fit postprocessing functions
get_y_rng <- function(fit, original_scale = TRUE) {
  yrng_done <- is_yrng_done(fit)
  if (!yrng_done) {
    stop("y rng was not done! you need options = list(do_yrng = TRUE)")
  }
  obs_model <- get_obs_model(fit)
  par_name <- if (obs_model == "gaussian") "y_rng_cont" else "y_rng_disc"
  out <- get_draws(fit, pars = par_name)
  if (obs_model == "gaussian" && original_scale) {
    scl <- dollar(fit@model@var_scalings, "y")
    out <- call_fun(scl@fun_inv, out)
  }
  colnames(out) <- NULL
  return(out)
}

#' Visualize the distribution of the obtained parameter draws
#'
#' @description
#' \itemize{
#'   \item \code{plot_draws} visualizes the distribution of any set of
#'   model parameters (defaults to kernel hyperparameters and possible
#'   observation model parameters)
#'   \item \code{plot_beta} visualizes the distribution of the
#'   individual-specific disease effect magnitude parameter draws
#'   \item \code{plot_effect_times} visualizes the input warping function for
#'   different parameter draws
#'   \item \code{plot_warp} visualizes the input warping function for
#'   different draws of the warping steepness parameter
#' }
#' @inheritParams get_draws
#' @param type plot type, allowed options are "intervals", "dens",
#' "areas", and "trace"
#' @param regex_pars regex for parameter names to plot
#' @param ... additional arguments for the \code{bayesplot} function
#' \code{\link[bayesplot]{mcmc_intervals}}, \code{\link[bayesplot]{mcmc_dens}},
#' \code{\link[bayesplot]{mcmc_areas}} or \code{\link[bayesplot]{mcmc_trace}}
#' @return a \code{ggplot} object or list of them
#' @name plot_draws
#' @family model fit visualization functions
NULL

#' @export
#' @rdname plot_draws
plot_draws <- function(fit,
                       type = "intervals",
                       regex_pars = c(
                         "alpha", "ell", "wrp",
                         "sigma", "phi", "gamma"
                       ),
                       ...) {
  check_type(fit, "lgpfit")
  allowed <- c("intervals", "areas", "trace", "dens")
  check_allowed(type, allowed)
  sf <- fit@stan_fit
  if (type == "dens") {
    h <- bayesplot::mcmc_dens(sf, regex_pars = regex_pars, ...)
  } else if (type == "trace") {
    h <- bayesplot::mcmc_trace(sf, regex_pars = regex_pars, ...)
  } else if (type == "areas") {
    h <- bayesplot::mcmc_areas(sf, regex_pars = regex_pars, ...)
  } else {
    h <- bayesplot::mcmc_intervals(sf, regex_pars = regex_pars, ...)
  }
  return(h)
}

#' @export
#' @rdname plot_draws
#' @inheritParams plot_inputwarp
#' @param num_points number of plot points
#' @param window_size width of time window
#' @param color_scheme deprecated argument, has no effect
#' @return a ggplot object
plot_warp <- function(fit, num_points = 300, window_size = 48,
                      color = colorset("red", "dark"), alpha = 0.5,
                      color_scheme = "brightblue") {
  check_type(fit, "lgpfit")
  R <- window_size
  num_ns <- dollar(fit@model@stan_input, "num_ns")
  dis_age <- seq(-R / 2, R / 2, length.out = num_points)
  out <- list()
  for (j in seq_len(num_ns)) {
    par_name <- paste0("wrp[", j, "]")
    draws <- rstan::extract(fit@stan_fit, pars = c(par_name))
    draws <- dollar(draws, par_name)
    out[[j]] <- plot_inputwarp(draws, dis_age, color, alpha)
  }

  # Return ggplot object or list of them
  L <- length(out)
  if (L == 0) stop("the model does not have warping parameters")
  simplify_list(out)
}

#' @export
#' @rdname plot_draws
plot_beta <- function(fit, type = "dens", ...) {
  check_type(fit, "lgpfit")
  num_heter <- fit@model@stan_input$num_heter
  if (num_heter == 0) {
    stop("there are no heterogeneous effects in the model")
  }
  h <- plot_draws(fit, type, regex_pars = "beta", ...)
  ptitle <- paste0(
    "Distribution of individual-specific effect magnitudes"
  )
  h <- h + ggplot2::ggtitle(label = ptitle)
  return(h)
}

#' @export
#' @rdname plot_draws
plot_effect_times <- function(fit, type = "areas", ...) {
  h <- plot_draws(fit, type, regex_pars = "teff[[]", ...)
  ptitle <- "Distribution of the inferred effect times"
  h <- h + ggplot2::ggtitle(label = ptitle)
  return(h)
}


#' Extract posterior draws
#'
#' @description These functions use \code{\link[rstan]{extract}} with
#' \code{permuted = FALSE} and \code{inc_warmup = FALSE}. Chains are stacked
#' so that the return value is always a 2-dimensional array.
#' \itemize{
#'  \item \code{get_draws} extracts posterior draws after warmup
#'   \item \code{get_draws.catch} function calls \code{get_draws} but catches
#' errors and returns NULL if an error occurs.
#'   \item \code{get_num_draws} returns totals number of post-warmup draws
#' }
#' @param fit an object of class \linkS4class{lgpfit}
#' @param draws Indices of parameter draws to return use. All post-warmup
#' draws are returned if this is \code{NULL}.
#' @param reduce Function used to reduce all parameter draws into
#' one set of parameters. Ignored if \code{NULL}, or if \code{draws} is not
#' \code{NULL}.
#' @param ... additional keyword arguments to \code{rstan::extract}
#' @return an array of shape \code{num_param_sets} x \code{num_params}
#' @family fit postprocessing functions
#' @name get_draws
NULL

#' @export
#' @rdname get_draws
get_draws <- function(fit, draws = NULL, reduce = NULL, ...) {
  check_type(fit, "lgpfit")
  s <- rstan::extract(fit@stan_fit, permuted = FALSE, inc_warmup = FALSE, ...)
  param_names <- dimnames(s)[[3]]
  s <- squeeze_second_dim(s) # squeeze the 'chains' dimension
  if (!is.null(draws)) {
    s <- s[draws, , drop = FALSE]
  } else {
    s <- apply_reduce(s, reduce)
  }
  colnames(s) <- param_names
  return(s)
}

#' @rdname get_draws
get_draws.catch <- function(fit, draws = NULL, reduce = NULL, ...) {
  tryCatch(
    {
      get_draws(fit, draws, reduce, ...)
    },
    error = function(e) {
      NULL
    }
  )
}

#' @export
#' @rdname get_draws
get_num_draws <- function(fit) {
  draws <- get_draws(fit, draws = NULL, reduce = NULL, pars = "alpha")
  nrow(draws)
}


#' Posterior summary
#'
#' @export
#' @inheritParams get_draws
#' @param ignore_pars names of parameters and generated quantities to ingore
#' @return a character representation
#' @family fit postprocessing functions
fit_summary <- function(fit,
                        ignore_pars = c(
                          "f_post", "f_latent", "eta",
                          "y_rng_disc", "y_rng_cont", "teff_raw", "lp__"
                        )) {
  check_type(fit, "lgpfit")
  print(fit@stan_fit, pars = ignore_pars, include = FALSE)
}
