#' Visualize a model posterior
#'
#' @param x an object of class \linkS4class{lgpfit}
#' @param y not used
#' @param ... keyword arguments passed to \code{\link{plot_posterior}}
#' @return a \code{ggplot} object
#' @family model fit visualization functions
setMethod(
  f = "plot",
  signature = signature(x = "lgpfit", y = "missing"),
  definition = function(x, ...) {
    plot_posterior(x, ...)
  }
)

#' Visualize posterior distribution of sampled parameters
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @param type plot type
#' @param regex_pars regex for parameter names to plot
#' @param ... other arguments to bayesplot functions
#' @return a \code{ggplot} object
#' @family model fit visualization functions
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
#' @family model fit visualization functions
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
#' @family model fit visualization functions
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
