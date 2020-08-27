#' Helper function
#'
#' @description Helper function for generic functions that work on
#' both of \linkS4class{lgpmodel} and \linkS4class{lgpfit} class objects.
#' @param object an object of class \linkS4class{lgpmodel} or
#' \linkS4class{lgpfit}
#' @return an object of class \linkS4class{lgpmodel}
object_to_model <- function(object) {
  allowed <- c("lgpmodel", "lgpfit")
  check_type(object, allowed)
  if (typeof(object) == "lgpfit") {
    out <- object@model
  } else {
    out <- object
  }
  return(out)
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
                           regex_pars = c("alpha", "ell", "wrp", "sigma", "phi"),
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

#' Visualize the input warping function for different parameter samples
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @param p number of plot points
#' @param R width of time window
#' @param color_scheme name of \code{bayesplot} color scheme
#' @return a \code{ggplot} object or list of them
plot_posterior_warp <- function(fit, p = 300, R = 48,
                                color_scheme = "brightblue") {
  check_type(fit, "lgpfit")

  # Colors
  scheme <- bayesplot::color_scheme_get(color_scheme)
  color_line <- scheme$dark
  color_inner <- scheme$light_highlight
  color_outer <- scheme$light

  # Plot
  num_ns <- fit@model@stan_input$num_ns
  ttt <- seq(-R / 2, R / 2, length.out = p)
  out <- list()
  for (j in seq_len(num_ns)) {
    par_name <- paste0("wrp[", j, "]")
    tsmr <- rstan::summary(fit@stan_fit, pars = c(par_name))$summary
    w_50 <- warp_input(ttt, a = tsmr[6])
    w_75 <- warp_input(ttt, a = tsmr[7])
    w_25 <- warp_input(ttt, a = tsmr[5])
    w_025 <- warp_input(ttt, a = tsmr[4])
    w_975 <- warp_input(ttt, a = tsmr[8])

    diseaseAge <- ttt
    DF <- data.frame(cbind(diseaseAge, w_50, w_75, w_25, w_025, w_975))

    # Create ggplot object
    h <- ggplot2::ggplot(DF, ggplot2::aes_string(x = "ttt", y = "w_50")) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "w_025", ymax = "w_975"),
        fill = color_outer
      ) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "w_25", ymax = "w_75"),
        fill = color_inner
      ) +
      ggplot2::geom_line(color = color_line)

    h <- h + ggplot2::labs(x = "Input", y = "Warped input")
    subt <- paste("Median steepness =", round(tsmr[6], 3))
    h <- h + ggplot2::ggtitle("Input-warping function", subtitle = subt)
    out[[j]] <- h
  }

  # Return ggplot object or list of them
  L <- length(out)
  if (L == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
}
