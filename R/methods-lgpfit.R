#' @describeIn lgpfit Print information and summary about the fit object.
setMethod("show", "lgpfit", function(object) {
  msg <- class_info("lgpfit")
  cat(msg)
  cat("\n")
  fit_summary(object)
})

#' @describeIn lgpfit Get the stored \linkS4class{lgpmodel} object.
#' Various properties of the returned object can be accessed as explained
#' in the documentation of \linkS4class{lgpmodel}.
setMethod("get_model", "lgpfit", function(object) {
  object@model
})

#' @describeIn lgpfit Get the stored \code{\link[rstan]{stanfit}} object.
#' Various properties of the returned object can be accessed or plotted
#' as explained
#' \href{https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html}{here}
#' or in the documentation of \code{\link[rstan]{stanfit}}.
setMethod("get_stanfit", "lgpfit", function(object) {
  object@stan_fit
})

#' @describeIn lgpfit Extract parameter draws. Uses \code{\link[rstan]{extract}} with
#' \code{permuted = FALSE} and \code{inc_warmup = FALSE}, so that the return
#' value is always a 2-dimensional array of shape
#' \code{num_param_sets} x \code{num_params}. Optional arguments
#' (\code{...}) are passed to \code{\link[rstan]{extract}}.
#' @param draws Indices of the parameter draws. \code{NULL} corresponds to
#' all post-warmup draws.
#' @param reduce Function used to reduce all parameter draws into
#' one set of parameters. Ignored if \code{NULL}, or if \code{draws} is not
#' \code{NULL}.
setMethod(
  "get_draws", "lgpfit",
  function(object, draws = NULL, reduce = NULL, ...) {
    s <- rstan::extract(object@stan_fit,
      permuted = FALSE,
      inc_warmup = FALSE,
      ...
    )
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
)

#' @describeIn lgpfit Visualize parameter draws. Optional arguments
#' (\code{...}) are passed to  \code{plotfun}.
#' @param x an \linkS4class{lgpfit} object to visualize
#' @param y unused argument
#' @param plotfun plotting function to use
#' @seealso For different plotting functions, see \code{\link{plot_draws}},
#' \code{\link{plot_beta}}, \code{\link{plot_warp}},
#' \code{\link{plot_effect_times}}
setMethod(
  "plot",
  signature = c("lgpfit", "missing"),
  function(x, y, plotfun = plot_draws, ...) {
    check_type(plotfun, "function")
    plotfun(x, ...)
  }
)

#' Print a fit summary.
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @param ignore_pars parameters and generated quantities to ignore from output
#' @returns \code{object} invisibly.
fit_summary <- function(fit,
                        ignore_pars = c(
                          "f_latent", "eta",
                          "teff_raw", "lp__"
                        )) {
  check_type(fit, "lgpfit")
  print(fit@stan_fit, pars = ignore_pars, include = FALSE)
}

#' Visualize the distribution of parameter draws
#'
#' @param fit an object of class \linkS4class{lgpfit}
#' @param type plot type, allowed options are "intervals", "dens",
#' "areas", and "trace"
#' @param regex_pars regex for parameter names to plot
#' @param ... additional arguments for the \code{bayesplot} function
#' \code{\link[bayesplot]{mcmc_intervals}}, \code{\link[bayesplot]{mcmc_dens}},
#' \code{\link[bayesplot]{mcmc_areas}} or \code{\link[bayesplot]{mcmc_trace}}
#' @return a \code{ggplot} object or list of them
#' @name plot_draws
NULL

#' @export
#' @describeIn plot_draws visualizes the distribution of any set of
#'   model parameters (defaults to kernel hyperparameters and possible
#'   observation model parameters)
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
#' @describeIn plot_draws visualizes the distribution of the
#'   individual-specific disease effect magnitude parameter draws
plot_beta <- function(fit, type = "dens", ...) {
  check_type(fit, "lgpfit")
  num_heter <- dollar(get_stan_input(fit), "num_heter")
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
#' @describeIn plot_draws visualizes the input warping function for
#'   different draws of the warping steepness parameter
#' @param num_points number of plot points
#' @param window_size width of time window
#' @param color line color
#' @param alpha line alpha
plot_warp <- function(fit, num_points = 300, window_size = 48,
                      color = colorset("red", "dark"), alpha = 0.5) {
  check_type(fit, "lgpfit")
  R <- window_size
  num_ns <- dollar(get_stan_input(fit), "num_ns")
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
  if (L == 0) stop("the model does not have input warping parameters")
  simplify_list(out)
}

#' @export
#' @describeIn plot_draws visualizes the input warping function for
#'   different parameter draws
plot_effect_times <- function(fit, type = "areas", ...) {
  num_uncrt <- dollar(get_stan_input(fit), "num_uncrt")
  if (num_uncrt == 0) {
    stop("there are no uncertain effect times in the model")
  }
  h <- plot_draws(fit, type, regex_pars = "teff[[]", ...)
  ptitle <- "Distribution of the inferred effect times"
  h <- h + ggplot2::ggtitle(label = ptitle)
  return(h)
}


determine_num_paramsets <- function(fit, draws, reduce) {
  # Decide number of output param sets based on total number of
  # posterior draws, possible reduction and subset of draw indices
  S <- fit@num_draws
  if (!is.null(reduce)) S <- 1
  if (!is.null(draws)) S <- length(draws)
  return(S)
}

#' Graphical posterior or prior predictive checks
#'
#' @export
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
  stop("not implemented")
  # y_rep <- get_y_rng(fit, original_scale = TRUE)
  # bayesplot::pp_check(y, y_rep, fun, ...)
}
