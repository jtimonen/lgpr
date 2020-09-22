#' Plot longitudinal data and/or model fit so that each subject/group has
#' their own panel
#'
#' @export
#' @description Data frames specified in arguments \code{df_data},
#' \code{df_signal}, \code{df_fit}, and \code{df_fit_err} must have a format
#' where
#' \itemize{
#'   \item the first column is the grouping factor (usually id)
#'   \item the second column is the x-axis variable (usually age)
#'   \item a column named \code{y} must contain the y-axis variable
#'   (not for \code{df_fit_err})
#'   \item a column named \code{lower} (\code{upper}) must contain the lower
#'   (upper) bound of error bar (only for \code{df_fit_err})
#'   \item a column named \code{draw} must be a factor that
#'   specifies the posterior draw using which the fit has been computed
#'   (only for \code{df_fit})
#' }
#' @param df_data A data frame containing the observations.
#' @param df_signal A data frame containing the true signal. Omitted if
#' \code{NULL}.
#' @param df_fit A data frame containing the model fit, or a list of data
#' frames. The list version can be used for example so that each list element
#' corresponds to the fit computed using one parameter draw. Omitted if
#' \code{NULL}.
#' @param df_fit_err A data frame containing error bars. Omitted if \code{NULL}.
#' Must be \code{NULL} if \code{df_fit} is a list.
#' @param teff_signal A named vector containing true effect times used to
#' generate the signal. Omitted if \code{NULL}.
#' @param teff_obs A named vector containing observed effect times. Omitted if
#' \code{NULL}.
#' @param i_test Indices of test points.
#' @param signal_color Line color for true signal.
#' @param fit_color Line color for model fit.
#' @param fit_err_color Color of the error ribbon.
#' @param vline_colors Two line colors for vertical lines
#' (true and obs. effect time).
#' @param fit_alpha Line opacity for model fit.
#' @param fit_err_alpha Opacity of the error ribbon.
#' @param nrow number of rows, an argument for
#' \code{\link[ggplot2]{facet_wrap}}
#' @param ncol number of columns, an argument for
#' \code{\link[ggplot2]{facet_wrap}}
#' @param y_transform A function.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @family plot APIs
plot_api_g <- function(df_data,
                       df_signal = NULL,
                       df_fit = NULL,
                       df_fit_err = NULL,
                       teff_signal = NULL,
                       teff_obs = NULL,
                       i_test = NULL,
                       signal_color = color_palette(2)[1],
                       fit_color = color_palette(2)[2],
                       fit_err_color = colorset("red", "light_highlight"),
                       vline_colors = colorset("gray", "mid_highlight"),
                       fit_alpha = 1.0,
                       fit_err_alpha = 1.0,
                       nrow = NULL,
                       ncol = NULL,
                       y_transform = function(x) x) {
  check_type(df_data, "data.frame")
  check_type(y_transform, "function")

  # Create the plot
  h <- plot_api_g_create(df_data, i_test, y_transform, nrow, ncol)
  group_by <- colnames(df_data)[1]
  x_name <- colnames(df_data)[2]

  # Add true and observed effect times
  col <- vline_colors
  h <- plot_api_g_add_effect_times(h, teff_signal, group_by, col, 1, 0.7)
  h <- plot_api_g_add_effect_times(h, teff_obs, group_by, col, 2, 0.7)

  # Add other layers
  h <- add_g_layer_ribbon(
    h, df_fit_err, x_name, group_by,
    fit_err_color, fit_err_alpha
  )
  h <- add_g_layer_fit(h, df_fit, x_name, group_by, fit_color, fit_alpha)
  h <- add_g_layer_signal(h, df_signal, x_name, group_by, signal_color)

  # Add data and return
  h <- h + ggplot2::geom_point()
  return(h)
}

#' Plot a generated/fit model component
#'
#' @export
#' @description Data frames specified in arguments \code{df},
#' and \code{df_err} must have a format where
#' \itemize{
#'   \item The first column is the grouping factor (usually id).
#'   \item The second column is the x-axis variable (usually age).
#'   \item The third column is the coloring factor. If name of the third
#'   column is \code{NA}, coloring is not done.
#'   \item A column named \code{y} must contain the y-axis variable
#'   (not for \code{df_err}).
#'   \item A column named \code{lower} (\code{upper}) must contain the lower
#'   (upper) bound of error bar (only for \code{df_err}).
#'   \item The posterior draw using which the fit has been computed can be
#'   specified with a factor named \code{_draw_} (only for \code{df}).
#' }
#' @param df a data frame
#' @param df_err a data frame
#' @param color_line line color
#' @param color_err ribbon color
#' @param alpha_line line opacity
#' @param alpha_err ribbon opacity
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @family plot APIs
plot_api_c <- function(df,
                       df_err = NULL,
                       color_line = color_palette(2)[2],
                       color_err = colorset("red", "light_highlight"),
                       alpha_line = 1.0,
                       alpha_err = 1.0) {

  # Create the plot
  group_by <- colnames(df)[1]
  x_name <- colnames(df)[2]
  color_by <- colnames(df)[3]
  color_by <- if (is.na(color_by)) NULL else color_by
  aes <- ggplot2::aes_string(x = x_name, y = "y")
  h <- ggplot2::ggplot(data = df, mapping = aes)

  # Add layers
  h <- add_c_layer_ribbon(h, df_err, x_name, group_by, color_by, alpha_err)
  h <- add_c_layer_line(h, df, x_name, group_by, color_by, alpha_line)
  return(h)
}

#' Initialize a grouped plot
#'
#' @inheritParams plot_api_g
#' @return a \code{ggplot} object
#' @family plot_api_g functions
plot_api_g_create <- function(df_data, i_test, y_transform, nrow, ncol) {

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
  aes <- plot_api_g_create_aes(df_data, type_name)
  h <- ggplot2::ggplot(df_data, aes)

  # Create faceting
  group_by <- colnames(df_data)[1]
  h <- add_g_layer_faceting(h, group_by, nrow, ncol)

  # Add shape and color scale
  if (show_test) {
    h <- h + ggplot2::scale_shape_manual(values = c(16, 16))
    col3 <- color_palette(3)[3]
    cols <- c(col3, "black")
    h <- h + ggplot2::scale_color_manual(values = cols)
  }
  return(h)
}

#' Helper function for plot_grouped
#'
#' @param data a data frame with three columns
#' @param type_name name of a factor that is to be distinguished by
#' plot marker and/or type
#' @return a ggplot aes object
#' @family plot_api_g functions
plot_api_g_create_aes <- function(data, type_name) {
  group_by <- colnames(data)[1]
  x_name <- colnames(data)[2]
  if (!is.null(type_name)) {
    aes <- ggplot2::aes_string(
      x = x_name,
      y = "y",
      group = group_by,
      pch = type_name,
      color = type_name
    )
  } else {
    aes <- ggplot2::aes_string(x = x_name, y = "y", group = group_by)
  }
  return(aes)
}

#' Add real or observed disease effect times to a grouped plot
#'
#' @param h a \code{ggplot} object
#' @param teff a named vector
#' @param group_by name of grouping variable
#' @param linecolor line color
#' @param linetype line type
#' @param lwd line width
#' @return updated \code{ggplot} object
#' @family plot_api_g functions
plot_api_g_add_effect_times <- function(h, teff, group_by,
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


#' Add layers to a grouped or componentwise plot
#'
#' @param h a ggplot object
#' @param df a data frame containing the things to add
#' @param x_name name of the x-axis variable
#' @param group_by name of the grouping factor
#' @param color_by name of the coloring factor
#' @param color color
#' @param alpha opacity
#' @return a modified ggplot object
#' @name add_layer
NULL

#' @rdname add_layer
add_c_layer_ribbon <- function(h, df, x_name, group_by, color_by, alpha) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    aes_rib <- ggplot2::aes_string(
      x = x_name,
      ymin = "lower",
      ymax = "upper",
      group = group_by,
      color = color_by,
      fill = color_by
    )
    h <- h + ggplot2::geom_ribbon(
      data = df,
      mapping = aes_rib,
      inherit.aes = FALSE,
      alpha = alpha
    )
  }
  return(h)
}

#' @rdname add_layer
add_c_layer_line <- function(h, df, x_name, group_by, color_by, alpha) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    if ("_draw_" %in% colnames(df)) {
      df <- add_factor_crossing(df, group_by, "_draw_", "group_x_draw")
      group_by <- "group_x_draw"
    }
    aes_fit <- ggplot2::aes_string(
      x = x_name,
      y = "y",
      color = color_by,
      group = group_by
    )
    h <- h + ggplot2::geom_line(
      data = df,
      mapping = aes_fit,
      inherit.aes = FALSE,
      alpha = alpha
    )
  }
  return(h)
}

#' @rdname add_layer
add_g_layer_ribbon <- function(h, df, x_name, group_by, color, alpha) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    aes_rib <- ggplot2::aes_string(
      x = x_name,
      ymin = "lower",
      ymax = "upper",
      group = group_by
    )
    h <- h + ggplot2::geom_ribbon(
      data = df,
      mapping = aes_rib,
      inherit.aes = FALSE,
      color = color,
      fill = color,
      alpha = alpha
    )
  }
  return(h)
}

#' @rdname add_layer
add_g_layer_fit <- function(h, df, x_name, group_by, color, alpha) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    if ("_draw_" %in% colnames(df)) {
      df <- add_factor_crossing(df, group_by, "_draw_", "group_x_draw")
      group_by <- "group_x_draw"
    }
    aes_fit <- ggplot2::aes_string(
      x = x_name,
      y = "y",
      group = group_by
    )
    h <- h + ggplot2::geom_line(
      data = df,
      mapping = aes_fit,
      inherit.aes = FALSE,
      color = color,
      alpha = alpha
    )
  }
  return(h)
}


#' @rdname add_layer
add_g_layer_signal <- function(h, df, x_name, group_by, color) {
  if (!is.null(df)) {
    check_type(df, "data.frame")
    aes_sig <- ggplot2::aes_string(
      x = x_name,
      y = "y",
      group = group_by
    )
    h <- h + ggplot2::geom_line(
      data = df,
      mapping = aes_sig,
      inherit.aes = FALSE,
      color = color
    )
  }
  return(h)
}

#' @rdname add_layer
#' @param facet_by factor to use for faceting
#' @inheritParams plot_api_g
add_g_layer_faceting <- function(h, facet_by, nrow, ncol) {
  f <- stats::as.formula(paste("~", facet_by))
  labfun <- ggplot2::label_both
  h <- h + ggplot2::facet_wrap(f, nrow = nrow, ncol = ncol, labeller = labfun)
  return(h)
}
