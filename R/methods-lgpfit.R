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
#' @param ... keyword arguments to \code{\link{plot_panel}}
#' @return a \code{ggplot object}
plot_fit <- function(fit, data, x_name = "age", y_name = "y",
                     group_by = "id", ...) {
  df_points <- data[c(group_by, x_name, y_name)]
  df_lines <- data[c(group_by, x_name)]
  f_post <- get_posterior_f(fit)
  lines <- list(
    mean = f_post$means$total,
    var = f_post$variances$total
  )

  h <- plot_panel(
    data = df_points,
    fit = list(x = df_lines, y = lines),
    ...
  )
  num_draws <- dim(lines$mean)[1]
  info <- paste("Showing analytic GP mean for", num_draws, "posterior draws.")
  h <- h + ggplot2::ggtitle("Model fit", subtitle = info)
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
#' @return a named list where each element is a named list
#' (length \code{num_components + 1}) of arrays of size
#' \code{num_draws} x \code{num_chains}
get_posterior_f <- function(fit) {
  check_type(fit, "lgpfit")
  f_sampled <- get_stan_input(fit)$is_f_sampled
  names <- get_component_names(fit)
  D <- length(names)
  R <- D + 1
  all_names <- c(names, "total")
  if (!f_sampled) {
    fp <- get_draws(fit, pars = "f_post")
    fp <- squeeze_second_dim(fp)
    alist <- array_to_arraylist(fp, 2 * R)
    m <- alist[1:R]
    v <- alist[(R + 1):(2 * R)]
    names(m) <- all_names
    names(v) <- all_names
    out <- list(means = m, variances = v)
  } else {
    fp <- get_draws(fit, pars = "f_latent")
    fp <- squeeze_second_dim(fp)
    alist <- array_to_arraylist(fp, R)
    names(alist) <- all_names
    out <- list(samples = alist)
  }
  return(out)
}
