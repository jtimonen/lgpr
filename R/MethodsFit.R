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

#' Extract posterior draws
#'
#' @export
#' @description Uses \code{rstan::extract} with \code{permuted = FALSE} and
#' \code{inc_warmup = FALSE}.
#' @param fit an object of class \linkS4class{lgpfit}
#' @param ... other keyword arguments to \code{rstan::extract}
#' @return a named list
#' @family fit postprocessing functions
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
#' @family fit postprocessing functions
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
#' @family fit postprocessing functions
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

#' Posterior summary
#'
#' @export
#' @param fit an object of class \linkS4class{lgpfit}
#' @param ignore_pars names of parameters and generated quantities to ingore
#' @return a character representation
#' @family fit postprocessing functions
fit_summary <- function(fit,
                        ignore_pars = c(
                          "f_post", "f_latent", "teff_raw",
                          "lp__"
                        )) {
  check_type(fit, "lgpfit")
  print(fit@stan_fit, pars = ignore_pars, include = FALSE)
}


#' Helper function
#'
#' @param fun a function
#' @param arr an array of shape \code{S} x \code{n}
#' @return an array with same shape as \code{arr}
scale_f_post_helper <- function(fun, arr) {
  check_type(fun, "function")
  DIM <- dim(arr)
  a <- as.numeric(t(arr))
  a <- fun(a)
  matrix(a, DIM[1], DIM[2], byrow = TRUE)
}

#' Helper function
#'
#' @param par_summary summary of warping parameter
#' @param dis_age the x-axis values
#' @param color_scheme name of \code{bayesplot} color scheme
#' @return a \code{ggplot} object
plot_posterior_warp_helper <- function(par_summary,
                                       dis_age,
                                       color_scheme = "brightblue") {

  # Get colors and quantiles
  scheme <- bayesplot::color_scheme_get(color_scheme)
  color_line <- scheme$dark
  color_inner <- scheme$light_highlight
  color_outer <- scheme$light
  w_50 <- warp_input(dis_age, a = par_summary[6])
  w_75 <- warp_input(dis_age, a = par_summary[7])
  w_25 <- warp_input(dis_age, a = par_summary[5])
  w_025 <- warp_input(dis_age, a = par_summary[4])
  w_975 <- warp_input(dis_age, a = par_summary[8])
  DF <- data.frame(cbind(dis_age, w_50, w_75, w_25, w_025, w_975))

  # Create ggplot object
  aes_line <- ggplot2::aes_string(x = "dis_age", y = "w_50")
  aes_outer <- ggplot2::aes_string(ymin = "w_025", ymax = "w_975")
  aes_inner <- ggplot2::aes_string(ymin = "w_25", ymax = "w_75")
  h <- ggplot2::ggplot(DF, aes_line) +
    ggplot2::geom_ribbon(aes_outer, fill = color_outer) +
    ggplot2::geom_ribbon(aes_inner, fill = color_inner) +
    ggplot2::geom_line(color = color_line)

  # Add titles and labels
  h <- h + ggplot2::labs(x = "Input", y = "Warped input")
  subt <- paste("Median steepness =", round(par_summary[6], 3))
  h <- h + ggplot2::ggtitle("Input-warping function", subtitle = subt)
  return(h)
}

#' Helper function
#'
#' @inheritParams plot_fit
#' @return a \code{ggplot object}
plot_fit_helper <- function(fit, data, x_name, y_name, group_by, draws) {
  df_data <- data[c(group_by, x_name, y_name)]
  df <- data[c(group_by, x_name)]
  f_post <- get_posterior_f(fit, draws)
  num_draws <- f_post$num_draws
  f_total <- f_post$f[["total"]]
  f_total <- scale_f_total(fit, f_total)

  # GP mean
  df_fit <- plot_fit_create_df(df, f_total$mean)

  # GP std
  if (num_draws > 1) {
    df_ribbon <- NULL
  } else {
    df_ribbon <- plot_fit_create_df_ribbon(df, f_total$mean, f_total$std, 2)
  }

  # Add effect time
  num_ns <- get_stan_input(fit)$num_ns
  if (num_ns == 0) {
    teff_obs <- NULL
  } else {
    da_name <- get_ns_covariates(fit)
    check_length(da_name, 1)
    teff_obs <- get_observed_effect_times(data, x_name, da_name, group_by)
  }

  # Plot title
  info <- "Showing analytic GP mean"
  if (num_draws > 1) {
    info <- paste(info, "for", num_draws, "posterior draws.")
  } else {
    info <- paste0(info, " and 2*std for posterior draw #", draws, ".")
  }

  # Return
  list(
    df_data = df_data,
    df_fit = df_fit,
    df_ribbon = df_ribbon,
    teff_obs = teff_obs,
    teff_fit = NULL,
    info = info,
    fit_alpha = line_alpha_fun(num_draws)
  )
}

#' Helper function
#'
#' @param df a data frame with group_by factor and x-variable
#' @param f an array of size \code{num_draws} x \code{num_obs}
#' @return a data frame
plot_fit_create_df <- function(df, f) {
  check_type(df, "data.frame")
  names <- colnames(df)
  S <- dim(f)[1]
  n <- dim(f)[2]
  X1 <- rep(df[, 1], S)
  X2 <- rep(df[, 2], S)
  X3 <- as.numeric(t(f))
  X4 <- rep(1:S, each = n)
  df_fit <- data.frame(as.factor(X1), X2, X3, as.factor(X4))
  colnames(df_fit) <- c(names, "f", "draw")
  return(df_fit)
}

#' Helper function
#'
#' @inheritParams plot_fit_create_df
#' @param f_mean an array of size \code{1} x \code{num_obs}
#' @param f_std an array of size \code{1} x \code{num_obs}
#' @param M multiplier for std
#' @return a data frame
plot_fit_create_df_ribbon <- function(df, f_mean, f_std, M) {
  check_type(df, "data.frame")
  check_not_null(M)
  f <- as.numeric(f_mean)
  std <- as.numeric(f_std)
  upper <- f + M * std
  lower <- f - M * std
  df_ribbon <- data.frame(lower, upper)
  df_ribbon <- cbind(df, df_ribbon)
  return(df_ribbon)
}
