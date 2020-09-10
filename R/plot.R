#' Visuzalize longitudinal data in panels
#'
#' @description Creates plots where each observation unit has a separate panel
#' @param df_data a data frame
#' @param df_signal a data frame or \code{NULL}
#' @param df_fit a data frame or \code{NULL}
#' @param df_ribbon a data frame or \code{NULL}
#' @param i_test indices of test points
#' @param linecolors two line colors (signal and fit)
#' @param vlinecolor color of vertical lines (true and obs. effect time)
#' @param ribbon_color ribbon color
#' @param nrow number of rows, an argument for \code{ggplot2::facet_wrap}
#' @param ncol number of columns, an argument for \code{ggplot2::facet_wrap}
#' @param y_transform a function
#' @param teff_true a named vector
#' @param teff_obs a named vector
#' @param teff_fit a list of named vectors
#' @param fit_alpha a value between 0 and 1
#' @param ribbon_alpha a value between 0 and 1
#' @return a \code{ggplot} object
plot_panel <- function(df_data,
                       df_signal = NULL,
                       df_fit = NULL,
                       df_ribbon = NULL,
                       i_test = NULL,
                       linecolors = color_palette(2),
                       vlinecolor = colorset("gray", "mid_highlight"),
                       ribbon_color = colorset("red", "light_highlight"),
                       nrow = NULL,
                       ncol = NULL,
                       y_transform = function(x) x,
                       teff_true = NULL,
                       teff_obs = NULL,
                       teff_fit = list(),
                       fit_alpha = 1.0,
                       ribbon_alpha = 1.0) {
  check_type(df_data, "data.frame")
  check_type(y_transform, "function")
  group_by <- colnames(df_data)[1]
  x_name <- colnames(df_data)[2]

  # Create the plot
  h <- plot_panel_create(df_data, i_test, y_transform, nrow, ncol)

  # Add true and observed effect times
  col <- vlinecolor
  h <- plot_panel_add_effect_times(h, teff_true, group_by, col, 1, 0.7)
  h <- plot_panel_add_effect_times(h, teff_obs, group_by, col, 2, 0.7)

  # Add ribbon
  if (!is.null(df_ribbon)) {
    check_type(df_ribbon, "data.frame")
    aes_rib <- ggplot2::aes_string(x = x_name, ymin = "lower", ymax = "upper")
    h <- h + ggplot2::geom_ribbon(
      data = df_ribbon,
      mapping = aes_rib,
      inherit.aes = FALSE,
      color = ribbon_color,
      fill = ribbon_color,
      alpha = ribbon_alpha
    )
  }

  # Add fit
  if (!is.null(df_fit)) {
    check_type(df_fit, "data.frame")
    aes_fit <- ggplot2::aes_string(x = x_name, y = "f", group = "draw")
    h <- h + ggplot2::geom_line(
      data = df_fit,
      mapping = aes_fit,
      inherit.aes = FALSE,
      color = linecolors[2],
      alpha = fit_alpha
    )
  }

  # Add signal
  if (!is.null(df_signal)) {
    check_type(df_signal, "data.frame")
    aes_sig <- plot_panel_create_aes(df_signal, NULL)
    h <- h + ggplot2::geom_line(
      data = df_signal,
      mapping = aes_sig,
      inherit.aes = FALSE,
      color = linecolors[1]
    )
  }

  # Add data
  h <- h + ggplot2::geom_point()

  return(h)
}

#' Initialize the panel plot
#'
#' @inheritParams plot_panel
#' @return a \code{ggplot} object
plot_panel_create <- function(df_data, i_test, y_transform, nrow, ncol) {

  # Create dummy variable for type (train/test)
  show_test <- !is.null(i_test)
  if (show_test) {
    n <- dim(df_data)[1]
    a <- rep("train", n)
    a[i_test] <- "test"
    type_name <- "Type"
    df_data[[type_name]] <- as.factor(a)
  } else {
    type_name <- NULL
  }

  # Create aes
  aes <- plot_panel_create_aes(df_data, type_name)
  h <- ggplot2::ggplot(df_data, aes)

  # Create faceting
  group_by <- colnames(df_data)[1]
  h <- plot_panel_add_faceting(h, group_by, nrow, ncol)

  # Add shape and color scale
  if (show_test) {
    h <- h + ggplot2::scale_shape_manual(values = c(16, 16))
    col3 <- color_palette(3)[3]
    cols <- c(col3, "black")
    h <- h + ggplot2::scale_color_manual(values = cols)
  }
  return(h)
}

#' Add faceting to the panel plot
#'
#' @inheritParams plot_panel_add_effect_times
#' @inheritParams plot_panel
#' @return a \code{ggplot} object
plot_panel_add_faceting <- function(h, group_by, nrow, ncol) {
  f <- stats::as.formula(paste("~", group_by))
  labfun <- ggplot2::label_both
  h <- h + ggplot2::facet_wrap(f, nrow = nrow, ncol = ncol, labeller = labfun)
  return(h)
}

#' Helper function
#'
#' @param data a data frame with three columns
#' @param type_name name of a factor that is to be distinguished by
#' plot marker and/or type
#' @return an aes object
plot_panel_create_aes <- function(data, type_name) {
  group_by <- colnames(data)[1]
  x_name <- colnames(data)[2]
  y_name <- colnames(data)[3]
  if (!is.null(type_name)) {
    aes <- ggplot2::aes_string(
      x = x_name,
      y = y_name,
      group = group_by,
      pch = type_name,
      color = type_name
    )
  } else {
    aes <- ggplot2::aes_string(x = x_name, y = y_name, group = group_by)
  }
  return(aes)
}

#' Add real or observed disease effect times to plot
#'
#' @param h a \code{ggplot} object
#' @param teff a named vector
#' @param group_by name of grouping variable
#' @param linecolor line color
#' @param linetype line type
#' @param lwd line width
#' @return updated \code{ggplot} object
plot_panel_add_effect_times <- function(h, teff, group_by,
                                        linecolor, linetype, lwd) {
  if (!is.null(teff)) {
    nams <- names(teff)
    z <- as.numeric(teff)
    df <- data.frame(z, nams)
    colnames(df) <- c("z", group_by)
    aes <- ggplot2::aes_string(xintercept = "z")
    h <- h + ggplot2::geom_vline(aes,
      na.rm = TRUE,
      data = df,
      color = linecolor,
      linetype = linetype,
      lwd = lwd
    )
  }
  return(h)
}


#' Plotting longitudinal data
#'
#' @export
#' @param data A data frame.
#' @param x_name Name of x-axis variable.
#' @param y_name Name of the y-axis variable.
#' @param group_by Name of grouping variable (must be a factor).
#' @param color_by Name of coloring variable (must be a factor).
#' @param facet_by Name of the faceting variable (must be a factor).
#' @param highlight Value of category of the \code{group_by} variable
#' that is highlighted. Can only be used if \code{color_by} is \code{NULL}.
#' @param main main plot title
#' @param sub plot subtitle
#' @return a \code{ggplot} object
plot_data <- function(data,
                      x_name = "age",
                      y_name = "y",
                      group_by = "id",
                      facet_by = NULL,
                      color_by = NULL,
                      highlight = NULL,
                      main = NULL,
                      sub = NULL) {
  check_type(data, "data.frame")

  # Create initial plot and add data
  df <- data[c(x_name, y_name, group_by, color_by, facet_by)]
  df <- plot_data_add_highlight_factor(df, group_by, highlight)
  skip_hl <- is.null(highlight)
  color_by <- if (skip_hl) color_by else paste0(group_by, "_")
  aes <- plot_data_aes(x_name, y_name, group_by, color_by)
  h <- ggplot2::ggplot(df, aes)

  # Add data
  h <- h + ggplot2::geom_line() + ggplot2::geom_point()

  # Add titles, faceting and coloring
  titles <- plot_data_titles(main, sub, data, group_by)
  label <- dollar(titles, "main")
  subtitle <- dollar(titles, "sub")
  h <- h + ggplot2::ggtitle(label = label, subtitle = subtitle)
  if (!is.null(facet_by)) {
    f <- stats::as.formula(paste("~", facet_by))
    h <- h + ggplot2::facet_wrap(f, labeller = ggplot2::label_both)
  }
  num_colors <- plot_data_num_colors(df, color_by)
  if (num_colors <= 4) {
    values <- color_palette(num_colors)
    h <- h + ggplot2::scale_color_manual(values = values)
  }

  return(h)
}

#' Create aes for plot_data
#'
#' @inheritParams plot_data
#' @return an \code{aes} object
plot_data_aes <- function(x_name, y_name, group_by, color_by) {
  if (is.null(color_by)) {
    aes <- ggplot2::aes_string(x = x_name, y = y_name, group = group_by)
  } else {
    aes <- ggplot2::aes_string(
      x = x_name, y = y_name, group = group_by,
      color = color_by
    )
  }
  return(aes)
}

#' Create titles for plot_data
#'
#' @inheritParams plot_data
#' @return a list
plot_data_titles <- function(main, sub, data, group_by) {
  g <- data[[group_by]]
  N <- length(unique(g))
  n <- dim(data)[1]
  if (is.null(main)) {
    main <- paste0(n, " data points")
  }
  if (is.null(sub)) {
    sub <- paste0(
      "Points that share the same value for '", group_by,
      "' are connected by a line (", N, " levels)."
    )
  }
  list(main = main, sub = sub, N = N, n = n)
}

#' Get number of distinct colors needsd by plot_data
#'
#' @inheritParams plot_data
#' @return a list
plot_data_num_colors <- function(data, color_by) {
  if (is.null(color_by)) {
    N <- 1
  } else {
    g <- data[[color_by]]
    N <- length(unique(g))
  }
  return(N)
}

#' Add factor to data frame for highlighting in plot
#'
#' @inheritParams plot_data
#' @param df data frame
#' @return a list
plot_data_add_highlight_factor <- function(df, group_by, highlight) {
  if (!is.null(highlight)) {
    check_length(highlight, 1)
    g <- df[[group_by]]
    s <- sum(g == highlight)
    if (s == 0) {
      str <- paste(levels(g), collapse = ", ")
      msg <- paste0(
        "Invalid <highlight> argument ", highlight, "! The ",
        "possible values of ", group_by, " are: {", str, "}."
      )
      stop(msg)
    }
    hl <- 1 + as.numeric(g == highlight)
    levels <- c("other", highlight)
    name <- paste0(group_by, "_")
    df[[name]] <- as.factor(levels[hl])
  }
  return(df)
}


#' Helper function
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
plot_fit_helper <- function(fit, draws) {
  df_data <- create_plot_df(fit)
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
    df_data = df_data,
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
