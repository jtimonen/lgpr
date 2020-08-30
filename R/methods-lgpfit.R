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
  df_data <- data[c(group_by, x_name, y_name)]
  df <- data[c(group_by, x_name)]
  list_f <- get_posterior_f(fit, draws)[["total"]]
  df_fit <- plot_fit_create_df(df, list_f)
  h <- plot_panel(df_data = df_data, df_fit = df_fit, ...)
  # TODO: edit
  num_draws <- dim(list_f$mean)[1]
  info <- paste("Showing analytic GP mean for", num_draws, "posterior draws.")
  h <- h + ggplot2::ggtitle("Model fit", subtitle = info)
  return(h)
}

#' Helper function
#'
#' @param df a data frame with group_by factor and x-variable
#' @param list_fit a list with fields mean and variance
#' @return a data frame
plot_fit_create_df <- function(df, list_fit) {
  check_type(df, "data.frame")
  check_type(list_fit, "list")
  names <- colnames(df)
  y_m <- list_fit$mean
  S <- dim(y_m)[1]
  n <- dim(y_m)[2]
  X1 <- rep(df[, 1], S)
  X2 <- rep(df[, 2], S)
  X3 <- as.numeric(t(y_m))
  X4 <- rep(1:S, each = n)
  df_fit <- data.frame(as.factor(X1), X2, X3, as.factor(X4))
  colnames(df_fit) <- c(names, "f", "draw")
  return(df_fit)
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
#' @param draws Indices of posterior draws for which to get \code{f}. This can
#' be a single integer, a vector of indices, or \code{NULL} (default). In the
#' latter case all draws are obtained.
#' @return Returns a named list of which has length equal to the number of
#' components plus one. Let \code{S = length(draws)}. Each list element is
#' \itemize{
#'   \item An array of size \code{S} x \code{num_obs}, if
#'   \code{is_sampled(model)} is \code{TRUE}. Each row of this array is one
#'   posterior draw of the function f.
#'   \item A list with fields \code{mean} and \code{variance}, if
#'   \code{is_sampled(model)} is \code{FALSE}. Both fields are arrays of size
#'   \code{S} x \code{num_obs}. These are the analytically computed means and
#'   variances for each posterior draw.
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
    mean <- alist[1:R]
    variance <- alist[(R + 1):(2 * R)]
    out <- zip_lists(mean, variance)
  } else {
    out <- array_to_arraylist(fp, R, draws)
  }
  names(out) <- all_names
  return(out)
}
