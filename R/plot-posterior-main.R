#' Deprecated function. Use plot_posterior_components instead.
#'
#' @export
#' @param fit fit
#' @param ... other arguments
#' @return Nothing. Throws an error when called.
plot_components_posterior <- function(fit, ...) {
  stop(
    "plot_components_posterior is deprecated. ",
    "Use plot_posterior_components instead."
  )
}

#' Plot posterior of f
#'
#' @export
#' @description This is a wrapper for \code{\link{plot_posterior_predictions}}.
#' @param fit An object of class \code{lgpfit}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param plot_uncertainty Should an uncertainty ribbon be plotted?
#' @param n_sds number of standard deviations for the uncertainty band width
#' @param data_marker pch for data points
#' @param ... additional arguments to \code{\link{plot_posterior_predictions}}
#'
#' @return a ggplot object
plot_posterior_f <- function(fit,
                             PRED = NULL,
                             plot_uncertainty = TRUE,
                             data_marker = 16,
                             n_sds = 2,
                             ...) {
  plot_obs_onset <- FALSE
  plot_t_effect_samples <- FALSE
  if (fit@model@stan_dat$UNCRT == 1) {
    plot_obs_onset <- TRUE
    plot_t_effect_samples <- TRUE
  }
  if (fit@model@info$sample_F) {
    alpha_line <- 0.1
  } else {
    alpha_line <- 1
  }
  h <- plot_posterior_predictions(fit,
    mode               = "posterior",
    PRED               = PRED,
    plot_uncertainty   = plot_uncertainty,
    n_sds              = n_sds,
    plot_obs_onset     = plot_obs_onset,
    plot_t_effect_samples = plot_t_effect_samples,
    alpha_line         = alpha_line,
    data_marker        = data_marker,
    ...
  )
  return(h)
}


#' Plot posterior predictive distribution
#'
#' @export
#' @description This is a wrapper for \code{\link{plot_posterior_predictions}}.
#' @param fit An object of class \code{lgpfit}.
#' @param PRED Predictions computed using \code{lgp_predict}.
#' @param uncertainty Either "none", "ribbon" or "errorbar".
#' @param test_data Test data set.
#' @param n_sds number of standard deviations for the uncertainty band width
#' @param data_marker pch for data points
#' @param ... additional arguments to \code{\link{plot_posterior_predictions}}
#' @return a ggplot object
plot_posterior_y <- function(fit, PRED,
                             uncertainty = "ribbon",
                             test_data = NULL,
                             data_marker = 16,
                             n_sds = 2,
                             ...) {
  if (uncertainty == "ribbon") {
    h <- plot_posterior_predictions(fit,
      mode = "predictive", PRED = PRED,
      test_data = test_data, n_sds = n_sds,
      data_marker = data_marker, ...
    )
  } else if (uncertainty == "errorbar") {
    h <- plot_posterior_predictions(fit,
      mode = "predictive", PRED = PRED,
      error_bar = TRUE, test_data = test_data,
      n_sds = n_sds, data_marker = data_marker, ...
    )
  } else {
    h <- plot_posterior_predictions(fit,
      mode = "predictive", PRED = PRED,
      plot_uncertainty = FALSE, test_data = test_data,
      n_sds = n_sds, data_marker = data_marker, ...
    )
  }
  return(h)
}


#' Visualize inferred components
#' @export
#' @inheritParams plot_posterior_components_sub1
#' @inheritParams plot_posterior_components_sub2
#' @return an object returned by \code{ggpubr::ggarrange} or a
#' list of ggplot2 objects
plot_posterior_components <- function(fit, subsamples = NULL,
                                      time_is_xvar = TRUE,
                                      PRED = NULL,
                                      marker = NULL,
                                      sample_idx = 1,
                                      n_sd = 2,
                                      ...) {
  if (is.null(PRED)) {
    if (is.null(marker) && is.null(subsamples)) {
      marker <- 16
    }
    h <- plot_posterior_components_sub1(
      fit, subsamples,
      time_is_xvar, marker, ...
    )
  } else {
    h <- plot_posterior_components_sub2(
      fit, PRED, sample_idx,
      time_is_xvar, n_sd, ...
    )
  }
  return(h)
}
