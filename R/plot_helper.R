#' Helper function for plot_warp
#'
#' @param par_summary summary of warping parameter
#' @param dis_age the x-axis values
#' @param color_scheme name of \code{bayesplot} color scheme
#' @return a \code{ggplot} object
plot_warp_helper <- function(par_summary,
                             dis_age,
                             color_scheme = "brightblue") {

  # Get colors and quantiles
  scheme <- bayesplot::color_scheme_get(color_scheme)
  color_line <- dollar(scheme, "dark")
  color_inner <- dollar(scheme, "light_highlight")
  color_outer <- dollar(scheme, "light")
  w_50 <- cpp_warp_input(dis_age, a = par_summary[6])
  w_75 <- cpp_warp_input(dis_age, a = par_summary[7])
  w_25 <- cpp_warp_input(dis_age, a = par_summary[5])
  w_025 <- cpp_warp_input(dis_age, a = par_summary[4])
  w_975 <- cpp_warp_input(dis_age, a = par_summary[8])
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
#' @param df_data a data frame returned by \code{\link{create_plot_df}}
#' @return a \code{ggplot object}
plot_fit_helper <- function(fit, df_data, draws) {
  df <- df_data[, 1:2]

  f_draws <- get_f(fit, draws)
  num_draws <- dollar(f_draws, "num_draws")
  f <- dollar(f_draws, "f")
  f_total <- dollar(f, "total")
  f_total <- scale_f_total(fit, f_total)

  # GP mean
  df_fit <- plot_fit_create_df(df, dollar(f_total, "mean"))

  # GP std
  if (num_draws > 1) {
    df_ribbon <- NULL
  } else {
    m <- dollar(f_total, "mean")
    s <- dollar(f_total, "std")
    df_ribbon <- plot_fit_create_df_ribbon(df, m, s, 2)
  }

  # Plot title
  info <- "Showing analytic GP mean"
  if (num_draws > 1) {
    info <- paste(info, "for", num_draws, " draws.")
  } else {
    info <- paste0(info, " and 2*std for draw #", draws, ".")
  }

  # Return
  list(
    df_fit = df_fit,
    df_ribbon = df_ribbon,
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
