#' Visualize posterior uncertainty in the disease effect times
#'
#' @export
#' @description Can only be used if the uncertainty of effect time was modeled.
#' @param fit An object of class \code{lgpfit}.
#' @param color_scheme Name of bayesplot color scheme.
#' @param prob Inner interval
#' @param prob_outer Outer interval
#' @param point_est Point estimate type
#' @family model fit visualization functions
#' @return a ggplot object
plot_effect_times <- function(fit,
                              color_scheme = "red",
                              prob = 1,
                              prob_outer = 1,
                              point_est = "none") {
  if (class(fit) != "lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  ptitle <- "Posterior distribution of the inferred disease effect time"
  sd <- fit@model@stan_dat
  if (sd$UNCRT == 0) {
    stop("The disease effect time was not modeled as uncertain!")
  }
  p <- plot_samples(fit,
    regex_pars = "T_effect",
    type = "areas",
    point_est = point_est,
    prob = prob,
    prob_outer = prob_outer,
    color_scheme = color_scheme
  )

  form <- fit@model@info$formula
  p <- p + ggplot2::labs(
    subtitle = paste("Model:", form),
    title = ptitle
  )
  return(p)
}


#' Visualize posterior samples of individual-specific disease effect magnitude
#' parameters
#'
#' @export
#' @description Can only be used if the disease effect was modeled
#' heterogeneously.
#' @param fit An object of class \code{lgpfit}.
#' @param color_scheme Name of bayesplot color scheme.
#' @param threshold Threshold for median.
#' @family model fit visualization functions
#' @return a ggplot object
plot_beta <- function(fit,
                      color_scheme = "red",
                      threshold = 0.5) {
  if (class(fit) != "lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  ptitle <- paste0(
    "Posterior distribution of individual-specific ",
    "disease effect magnitudes"
  )
  aff <- affected(fit, threshold = threshold)
  df <- as.data.frame(fit@stan_fit)
  ibeta <- grep("beta", names(df))
  df <- df[, ibeta]
  colnames(df) <- paste("id = ", names(aff), sep = "")
  bayesplot::color_scheme_set(scheme = color_scheme)
  p <- bayesplot::mcmc_dens(df)
  beta <- "beta" # avoid note from R CMD check
  p <- p + ggplot2::xlab(expression(beta))

  str <- paste("Affected individuals: ",
    paste(names(which(aff)), collapse = ", "),
    sep = ""
  )
  p <- p + ggplot2::labs(
    subtitle = str,
    title = ptitle
  )
  return(p)
}
