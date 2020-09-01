#' Visualize a model posterior
#'
#' @param x an object of class \linkS4class{lgpfit}
#' @param y not used
#' @param ... keyword arguments passed to \code{\link{plot_posterior}}
#' @return a \code{ggplot} object
setMethod(
  f = "plot",
  signature = signature(x = "lgpfit", y = "missing"),
  definition = function(x, ...) {
    plot_posterior(x, ...)
  }
)

#' Posterior summary
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @return a character representation
fit_summary <- function(fit) {
  check_type(fit, "lgpfit")
  print(fit@stan_fit, pars = c("f_post", "lp__"), include = FALSE)
}

#' Visualize posterior distribution of sampled parameters
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @param type plot type
#' @param regex_pars regex for parameter names to plot
#' @param ... other arguments to bayesplot functions
#' @return a \code{ggplot} object
plot_posterior <- function(fit,
                           type = "intervals",
                           regex_pars = c(
                             "alpha", "ell", "wrp",
                             "sigma", "phi"
                           ),
                           ...) {
  check_type(fit, "lgpfit")
  sf <- fit@stan_fit
  if (type == "dens") {
    h <- bayesplot::mcmc_dens(sf, regex_pars = regex_pars, ...)
  } else if (type == "trace") {
    h <- bayesplot::mcmc_trace(sf, regex_pars = regex_pars, ...)
  } else {
    h <- bayesplot::mcmc_intervals(sf, regex_pars = regex_pars, ...)
  }
  return(h)
}

#' Visualize a model fit against longitudinal data set
#'
#' @export
#' @description Creates plots where each observation unit has a separate panel.
#' @param fit an object of class \linkS4class{lgpfit}
#' @param data a data frame
#' @param x_name name of x-axis variable
#' @param y_name name of y-axis variable
#' @param group_by grouping variable
#' @param draws see the \code{draws} argument of \code{\link{get_posterior_f}}
#' @param ... keyword arguments to \code{\link{plot_panel}}
#' @return a \code{ggplot object}
plot_fit <- function(fit, data, x_name = "age", y_name = "y",
                     group_by = "id", draws = NULL, ...) {
  DF <- plot_fit_helper(fit, data, x_name, y_name, group_by, draws)
  h <- plot_panel(
    df_data = DF$df_data,
    df_fit = DF$df_fit,
    df_ribbon = DF$df_ribbon,
    fit_alpha = DF$fit_alpha,
    teff_obs = DF$teff_obs,
    teff_fit = DF$teff_fit,
    ...
  )
  h <- h + ggplot2::ggtitle("Model fit", subtitle = DF$info)
  return(h)
}

#' Visualize the input warping function for different parameter samples
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @param p number of plot points
#' @param R width of time window
#' @inheritParams plot_posterior_warp_helper
#' @return a \code{ggplot} object or list of them
plot_posterior_warp <- function(fit, p = 300, R = 48,
                                color_scheme = "brightblue") {
  check_type(fit, "lgpfit")
  num_ns <- fit@model@stan_input$num_ns
  dis_age <- seq(-R / 2, R / 2, length.out = p)
  out <- list()
  for (j in seq_len(num_ns)) {
    par_name <- paste0("wrp[", j, "]")
    par_summary <- rstan::summary(fit@stan_fit, pars = c(par_name))$summary
    out[[j]] <- plot_posterior_warp_helper(par_summary, dis_age, color_scheme)
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

#' Extract posterior draws
#'
#' @description Uses \code{rstan::extract} with \code{permuted = FALSE} and
#' \code{inc_warmup = FALSE}.
#' @param fit an object of class \linkS4class{lgpfit}
#' @param ... other keyword arguments to \code{rstan::extract}
#' @return a named list
get_draws <- function(fit, ...) {
  check_type(fit, "lgpfit")
  rstan::extract(fit@stan_fit, permuted = FALSE, inc_warmup = FALSE, ...)
}

#' Extract posterior of the function f and its components
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
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
get_posterior_f <- function(fit, draws = NULL) {
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
    f_out <- array_to_arraylist(fp, R, draws)
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
#' @param fit an object of class \linkS4class{lgpfit}
#' @param f_total a list with fields \code{mean} and \code{std}
#' @return a similar object as \code{f_total}
scale_f_total <- function(fit, f_total) {
  check_type(fit, "lgpfit")
  check_not_null(f_total)
  f_sampled <- is_f_sampled(fit)
  check_false(f_sampled)
  fun_inv <- fit@model@var_scalings$y@fun_inv
  f_total$mean <- scale_f_post_helper(fun_inv, f_total$mean)
  f_total$std <- scale_f_post_helper(fun_inv, f_total$std)
  return(f_total)
}
