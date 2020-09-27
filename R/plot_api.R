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
#' @param df A data frame containing the model fit, or a list of data
#' frames. The list version can be used for example so that each list element
#' corresponds to the fit computed using one parameter draw. Omitted if
#' \code{NULL}.
#' @param df_err A data frame containing error bars. Omitted if \code{NULL}.
#' Must be \code{NULL} if \code{df_fit} is a list.
#' @param teff_signal A named vector containing true effect times used to
#' generate the signal. Omitted if \code{NULL}.
#' @param teff_obs A named vector containing observed effect times. Omitted if
#' \code{NULL}.
#' @param i_test Indices of test points.
#' @param color_signal Line color for true signal.
#' @param color Line color for model fit.
#' @param color_err Color of the error ribbon.
#' @param color_vlines Two line colors for vertical lines
#' (true and obs. effect time).
#' @param alpha Line opacity for model fit.
#' @param alpha_err Opacity of the error ribbon.
#' @param nrow number of rows, an argument for
#' \code{\link[ggplot2]{facet_wrap}}
#' @param ncol number of columns, an argument for
#' \code{\link[ggplot2]{facet_wrap}}
#' @param y_transform A function.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @family plot APIs
plot_api_g <- function(df_data,
                       df_signal = NULL,
                       df = NULL,
                       df_err = NULL,
                       teff_signal = NULL,
                       teff_obs = NULL,
                       i_test = NULL,
                       color_signal = color_palette(2)[1],
                       color = color_palette(2)[2],
                       color_err = colorset("red", "light_highlight"),
                       color_vlines = colorset("gray", "mid_highlight"),
                       alpha = 1.0,
                       alpha_err = 0.5,
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
  col <- color_vlines
  h <- plot_api_g_add_effect_times(h, teff_signal, group_by, col, 1, 0.7)
  h <- plot_api_g_add_effect_times(h, teff_obs, group_by, col, 2, 0.7)

  # Add other layers
  h <- add_g_layer_ribbon(
    h, df_err, x_name, group_by,
    color_err, alpha_err
  )
  h <- add_g_layer_fit(h, df, x_name, group_by, color, alpha)
  h <- add_g_layer_signal(h, df_signal, x_name, group_by, color_signal)

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
#' @param alpha line opacity
#' @param alpha_err ribbon opacity
#' @param no_err hide error bar even when it would normally be plotted?
#' @param no_line hide line even when it would normally be plotted?
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @family plot APIs
plot_api_c <- function(df,
                       df_err = NULL,
                       alpha = 1.0,
                       alpha_err = 0.2,
                       no_err = FALSE,
                       no_line = FALSE) {

  # Create the plot
  group_by <- colnames(df)[1]
  x_name <- colnames(df)[2]
  color_by <- colnames(df)[3]
  color_by <- if (is.na(color_by)) NULL else color_by
  aes <- ggplot2::aes_string(x = x_name, y = "y")
  h <- ggplot2::ggplot(data = df, mapping = aes)

  # Add layers
  if (!no_err) {
    h <- add_c_layer_ribbon(h, df_err, x_name, group_by, color_by, alpha_err)
  }
  if (!no_line) {
    h <- add_c_layer_line(h, df, x_name, group_by, color_by, alpha)
  }
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
    if (!is.null(color_by)) {
      h <- h + scale_color(5)
    }
  }
  return(h)
}

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
      alpha = alpha,
      color = NA # no ribbon edge
    )
    if (!is.null(color_by)) {
      h <- h + scale_fill(5)
    }
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
      color = NA,
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
