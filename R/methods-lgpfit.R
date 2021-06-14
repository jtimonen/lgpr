#' @export
#' @describeIn lgpfit Print information and summary about the fit object.
setMethod("show", "lgpfit", function(object) {
  desc <- class_info("lgpfit")
  cat(desc)
  cat("\n")
  fit_summary(object)
})

#' @export
#' @describeIn lgpfit Get names of model components.
setMethod("component_names", "lgpfit", function(object) {
  component_names(get_model(object))
})

#' @export
#' @describeIn lgpfit Apply postprocessing. Returns an updated
#' \linkS4class{lgpfit} object (copies data).
#' @param verbose Can the method print any messages?
setMethod("postproc", "lgpfit", function(object, verbose = TRUE) {

  # Compute pred that can be used to compute relevances
  if (contains_postproc(object)) {
    msg <- paste0(
      "Object already contains postprocessing information!",
      " You can remove it by calling clear_postproc()."
    )
    stop(msg)
  }
  pred <- pred(fit = object, reduce = NULL, verbose = verbose)

  # Return
  new("lgpfit",
    model = object@model,
    stan_fit = object@stan_fit,
    num_draws = object@num_draws,
    postproc_results = list(pred = pred)
  )
})

#' @export
#' @describeIn lgpfit Check if object contains postprocessing information.
setMethod("contains_postproc", "lgpfit", function(object) {
  length(object@postproc_results) > 0
})

#' @export
#' @describeIn lgpfit Returns an updated (copies data)
#' \linkS4class{lgpfit} object without any postprocessing information.
setMethod("clear_postproc", "lgpfit", function(object) {
  new("lgpfit",
    model = object@model,
    stan_fit = object@stan_fit,
    num_draws = object@num_draws,
    postproc_results = list()
  )
})

#' @export
#' @describeIn lgpfit Get the stored \linkS4class{lgpmodel} object.
#' Various properties of the returned object can be accessed as explained
#' in the documentation of \linkS4class{lgpmodel}.
setMethod("get_model", "lgpfit", function(object) {
  object@model
})

#' @export
#' @describeIn lgpfit Get the stored \code{\link[rstan]{stanfit}} object.
#' Various properties of the returned object can be accessed or plotted
#' as explained
#' \href{https://CRAN.R-project.org/package=rstan/vignettes/stanfit-objects.html}{here}
#' or in the documentation of \code{\link[rstan]{stanfit}}.
setMethod("get_stanfit", "lgpfit", function(object) {
  object@stan_fit
})

#' @export
#' @describeIn lgpfit Determine if inference was done by sampling
#' the latent signal \code{f} (and its components).
setMethod("is_f_sampled", "lgpfit", function(object) {
  is_f_sampled(object@model)
})

#' @export
#' @describeIn lgpfit Visualize parameter draws using \code{\link{plot_draws}}.
#' @param x an \linkS4class{lgpfit} object to visualize
#' @param y unused argument
#' @seealso For more detailed plotting functions, see \code{\link{plot_draws}},
#' \code{\link{plot_beta}}, \code{\link{plot_warp}},
#' \code{\link{plot_effect_times}}
setMethod(
  "plot",
  signature = c("lgpfit", "missing"),
  function(x, y) {
    plot_draws(fit = x)
  }
)


#' Extract parameter draws from lgpfit or stanfit
#'
#' @export
#' @description Uses \code{\link[rstan]{extract}}
#' with \code{permuted = FALSE} and \code{inc_warmup = FALSE}.
#' @param object An object of class \linkS4class{lgpfit} or \code{stanfit}.
#' @param draws Indices of the parameter draws. \code{NULL} corresponds to
#' all post-warmup draws.
#' @param reduce Function used to reduce all parameter draws into
#' one set of parameters. Ignored if \code{NULL}, or if \code{draws} is not
#' \code{NULL}.
#' @param ... Additional arguments to \code{rstan::extract()}.
#' @return The return value is always a 2-dimensional array of shape
#' \code{num_param_sets} x \code{num_params}.
#' @family main functions
get_draws <- function(object, draws = NULL, reduce = NULL, ...) {
  if (!is(object, "stanfit")) {
    sf <- get_stanfit(object)
  } else {
    sf <- object
  }
  s <- rstan::extract(sf, permuted = FALSE, inc_warmup = FALSE, ...)
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
  sf <- get_stanfit(fit)
  print(sf, pars = ignore_pars, include = FALSE)
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
#' @param verbose Can any output be printed?
#' @name plot_draws
NULL

#' @export
#' @describeIn plot_draws visualizes the distribution of any set of
#'   model parameters (defaults to kernel hyperparameters and possible
#'   observation model parameters)
#' @family main plot functions
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
  sf <- get_stanfit(fit)
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
plot_beta <- function(fit, type = "dens", verbose = TRUE, ...) {
  check_type(fit, "lgpfit")
  num_heter <- dollar(get_stan_input(fit), "num_heter")
  if (num_heter == 0) {
    stop("there are no heterogeneous effects in the model")
  }
  if (verbose) {
    print(beta_teff_idx_info(fit)) # print parameter indices
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
  sf <- get_stanfit(fit)
  for (j in seq_len(num_ns)) {
    par_name <- paste0("wrp[", j, "]")
    draws <- rstan::extract(sf, pars = c(par_name))
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
plot_effect_times <- function(fit, type = "areas", verbose = TRUE, ...) {
  num_uncrt <- dollar(get_stan_input(fit), "num_uncrt")
  if (num_uncrt == 0) {
    stop("there are no uncertain effect times in the model")
  }
  if (verbose) {
    print(beta_teff_idx_info(fit)) # print parameter indices
  }
  h <- plot_draws(fit, type, regex_pars = "teff[[]", ...)
  ptitle <- "Distribution of the inferred effect times"
  h <- h + ggplot2::ggtitle(label = ptitle)
  return(h)
}

determine_num_paramsets <- function(stan_fit, draws, reduce) {
  # Decide number of output param sets based on total number of
  # posterior draws, possible reduction and subset of draw indices
  stopifnot(is(stan_fit, "stanfit"))
  S <- get_num_postwarmup_draws(stan_fit)
  if (!is.null(reduce)) S <- 1
  if (!is.null(draws)) S <- length(draws)
  return(S)
}
