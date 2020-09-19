#' Posterior or prior predictive checks for lgpfit objects
#'
#' @export
#' @inheritParams get_draws
#' @param data the original data frame
#' @param fun \code{bayesplot} function name
#' @param ... arguments passed to the default [bayesplot::pp_check()] method in
#' \code{bayesplot}
#' @return a \code{ggplot} object
ppc <- function(fit, data, fun = bayesplot::ppc_dens_overlay, ...) {
  check_type(fit, "lgpfit")
  check_type(data, "data.frame")
  check_type(fun, "function")
  y_name <- get_y_name(fit)
  y <- dollar(data, y_name)
  y_rep <- get_y_rng(fit, original_scale = TRUE)
  bayesplot::pp_check(y, y_rep, fun, ...)
}

#' Visualize a model fit against longitudinal data set
#'
#' @export
#' @description Creates plots where each observation unit has a separate panel.
#' @inheritParams get_draws
#' @inheritParams create_plot_df
#' @param draws see the \code{draws} argument of \code{\link{get_f}}
#' @param ... keyword arguments to \code{\link{plot_panel}}
#' @return a \code{ggplot} object
#' @family model fit visualization functions
plot_fit <- function(fit, x = "age", group_by = "id", draws = NULL, ...) {
  df_data <- create_plot_df(fit, x, group_by)
  DF <- plot_fit_helper(fit, df_data, draws)
  h <- plot_panel(
    df_data = df_data,
    df_fit =  dollar(DF, "df_fit"),
    df_ribbon = dollar(DF, "df_ribbon"),
    fit_alpha = dollar(DF, "fit_alpha"),
    ...
  )
  h <- h + ggplot2::ggtitle("Model fit", subtitle = DF$info)
  return(h)
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
#' @param num_points number of plot points
#' @param window_size width of time window
#' @param color_scheme name of \code{bayesplot} color scheme
plot_warp <- function(fit, num_points = 300, window_size = 48,
                      color_scheme = "brightblue") {
  check_type(fit, "lgpfit")
  R <- window_size
  num_ns <- dollar(fit@model@stan_input, "num_ns")
  dis_age <- seq(-R / 2, R / 2, length.out = num_points)
  out <- list()
  for (j in seq_len(num_ns)) {
    par_name <- paste0("wrp[", j, "]")
    summary <- rstan::summary(fit@stan_fit, pars = c(par_name))
    par_summary <- dollar(summary, "summary")
    out[[j]] <- plot_warp_helper(par_summary, dis_age, color_scheme)
  }

  # Return ggplot object or list of them
  L <- length(out)
  if (L == 1) {
    return(out[[1]])
  } else {
    if (L == 0) {
      stop("The model does not have warping parameters.")
    }
    return(out)
  }
}

#' @export
#' @rdname plot_draws
plot_beta <- function(fit, type = "dens", ...) {
  h <- plot_draws(fit, type, regex_pars = "beta", ...)
  ptitle <- paste0(
    "Distribution of individual-specific ",
    "disease effect magnitudes"
  )
  h <- h + ggplot2::ggtitle(label = ptitle)
  return(h)
}

#' @export
#' @rdname plot_draws
plot_effect_times <- function(fit, type = "areas", ...) {
  h <- plot_draws(fit, type, regex_pars = "teff[[]", ...)
  ptitle <- "Distribution of the inferred disease effect times"
  h <- h + ggplot2::ggtitle(label = ptitle)
  return(h)
}


#' Extract posterior draws
#'
#' @export
#' @description These functions use \code{\link[rstan]{extract}} with
#' \code{permuted = FALSE} and \code{inc_warmup = FALSE}.
#' The \code{get_draws_catch} function calls \code{get_draws} but catches
#' errors and returns NULL if an error occurs.
#' @param fit an object of class \linkS4class{lgpfit}
#' @param stack_chains should this return a 2-dimensional array?
#' @param ... additional keyword arguments to \code{rstan::extract}
#' @return a 3 or 2-dimensional array
#' @family fit postprocessing functions
#' @name get_draws
NULL

#' @export
#' @rdname get_draws
get_draws <- function(fit, stack_chains = FALSE, ...) {
  check_type(fit, "lgpfit")
  s <- rstan::extract(fit@stan_fit, permuted = FALSE, inc_warmup = FALSE, ...)
  s <- if (stack_chains) squeeze_second_dim(s) else s
  return(s)
}

#' @export
#' @rdname get_draws
get_draws_catch <- function(fit, stack_chains = FALSE, ...) {
  tryCatch(
    {
      get_draws(fit, stack_chains, ...)
    },
    error = function(e) {
      NULL
    }
  )
}

#' @export
#' @rdname get_draws
get_num_draws <- function(fit) {
  draws <- get_draws(fit, stack_chains = TRUE, pars = "alpha")
  nrow(draws)
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
  f_sampled <- is_f_sampled(fit)
  if (!f_sampled) stop("f was not sampled!")
  obs_model <- get_obs_model(fit)
  par_name <- if (obs_model == "gaussian") "y_rng_cont" else "y_rng_disc"
  out <- get_draws(fit, pars = par_name)
  out <- squeeze_second_dim(out)
  if (obs_model == "gaussian" && original_scale) {
    scl <- dollar(fit@model@var_scalings, "y")
    out <- call_fun(scl@fun_inv, out)
  }
  colnames(out) <- NULL
  return(out)
}

#' Extract draws of the function f and its components
#'
#' @export
#' @inheritParams get_draws
#' @param draws Indices of posterior draws for which to get \code{f}. This can
#' be a single integer, a vector of indices, or \code{NULL} (default). In the
#' latter case all draws are obtained.
#' @return Returns a list with names \code{num_draws} and \code{f}.
#' The latter is a named list of which has length equal to the number of
#' components plus one. Let \code{S = length(draws)}. Each list element is
#' \itemize{
#'   \item Array of size \code{S} x \code{num_obs}, where each row is one
#'   posterior draw of the function f, \code{is_sampled(model)} is \code{TRUE}.
#'   \item A list with fields \code{mean} and \code{std}, if
#'   \code{is_sampled(model)} is \code{FALSE}. Both fields are arrays of size
#'   \code{S} x \code{num_obs}. These are the analytically computed means and
#'   standard deviations for each posterior draw.
#' }
#' @family fit postprocessing functions
get_f <- function(fit, draws = NULL) {
  check_type(fit, "lgpfit")
  f_sampled <- is_f_sampled(fit)
  names <- get_component_names(fit)
  D <- length(names)
  R <- D + 1
  all_names <- c(names, "total")
  pars <- if (f_sampled) "f_latent" else "f_post"
  fp <- get_draws(fit, pars = pars)
  fp <- squeeze_second_dim(fp)
  S <- dim(fp)[1]
  if (is.null(draws)) {
    draws <- c(1:S)
  }
  if (!f_sampled) {
    alist <- array_to_arraylist(fp, 2 * R, draws)
    mean <- alist[1:R] # means
    std <- alist[(R + 1):(2 * R)] # stds
    f_out <- zip_lists(mean, std)
  } else {
    f_out <- array_to_arraylist(fp, D, draws)
    f_out <- add_sum_arraylist(f_out)
  }
  names(f_out) <- all_names

  # Return
  list(
    f = f_out,
    num_draws = length(draws)
  )
}

#' Scale the function f posterior to original unnormalized scale
#'
#' @export
#' @description Can only be used with Gaussian observation model.
#' @inheritParams get_draws
#' @param f_total a list with fields \code{mean} and \code{std}
#' @return a similar object as \code{f_total}
#' @family fit postprocessing functions
scale_f_total <- function(fit, f_total) {
  check_type(fit, "lgpfit")
  check_not_null(f_total)
  f_sampled <- is_f_sampled(fit)
  check_false(f_sampled)
  y_scl <- dollar(fit@model@var_scalings, "y")
  fun_inv <- y_scl@fun_inv
  f_total$mean <- scale_f_helper(fun_inv, f_total$mean)
  f_total$std <- scale_f_helper(fun_inv, f_total$std)
  return(f_total)
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

#' Helper function
#'
#' @param fun a function
#' @param arr an array of shape \code{S} x \code{n}
#' @return an array with same shape as \code{arr}
scale_f_helper <- function(fun, arr) {
  check_type(fun, "function")
  DIM <- dim(arr)
  a <- as.numeric(t(arr))
  a <- fun(a)
  matrix(a, DIM[1], DIM[2], byrow = TRUE)
}

#' Extract modeled signal or its components
#'
#' @description
#' \itemize{
#'   \item \code{get_f_total} gets total signal for each sample
#'   \item \code{get_f_components} gets signal components for each sample
#'   \item \code{get_g} gets total signal, mapped through
#'   the inverse link function, for each sample
#' }
#' @param field Should be \code{"mean"}, \code{"std"}, or \code{NULL}.
#' @inheritParams get_draws
#' @return an array of shape \code{num_samples} x \code{num_obs} (or a list of
#' them for \code{get_f_components})
#' @name signal
NULL

#' @rdname signal
get_f_total <- function(fit, field = "mean") {
  f <- dollar(get_f(fit), "f")
  f <- dollar(f, "total")
  if (is_f_sampled(fit)) {
    check_null(field, "f was sampled in this model")
  } else {
    allowed <- c("mean", "std")
    check_allowed(field, allowed)
    f <- dollar(f, field)
  }
  colnames(f) <- NULL
  return(f)
}

#' @rdname signal
get_f_components <- function(fit, field = "mean") {
  f <- dollar(get_f(fit), "f")
  f[["total"]] <- NULL
  if (is_f_sampled(fit)) {
    check_null(field, "f was sampled in this model")
  } else {
    allowed <- c("mean", "std")
    check_allowed(field, allowed)
    f <- lapply(f, "[[", field)
  }
  return(f)
}

#' @rdname signal
get_g <- function(fit) {
  flag <- is_f_sampled(fit)
  field <- if (flag) NULL else "mean"
  f <- get_f_total(fit, field)
  likelihood <- get_obs_model(fit)
  link_inv(f, likelihood)
}
