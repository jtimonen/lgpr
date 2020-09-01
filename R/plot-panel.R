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
    names <- names(teff)
    z <- as.numeric(teff)
    df <- data.frame(z, names)
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
