#' Visuzalize longitudinal data in panels
#'
#' @description Creates plots where each observation unit has a separate panel
#' @param data a data frame
#' @param signal a data frame
#' @param fit a list of data frames
#' @param i_test indices of test points
#' @param linecolors three line colors
#' @param nrow number of rows, an argument for \code{ggplot2::facet_wrap}
#' @param ncol number of columns, an argument for \code{ggplot2::facet_wrap}
#' @param y_transform a function
#' @param true_teff a named vector
#' @param signal_teff a named vector
#' @param fit_teff a list of named vectors
#' @return a \code{ggplot} object
plot_panel <- function(data,
                       signal = NULL,
                       fit = NULL,
                       i_test = NULL,
                       linecolors = color_palette(2),
                       nrow = NULL,
                       ncol = NULL,
                       y_transform = function(x) x,
                       true_teff = NULL,
                       signal_teff = NULL,
                       fit_teff = list()) {
  check_type(data, "data.frame")
  check_type(y_transform, "function")
  group_by <- colnames(data)[1]

  # Create the plot
  h <- plot_panel_create(data, i_test, y_transform, nrow, ncol)

  # Add true and observed effect times
  col1 <- colorset("red")
  col2 <- colorset("red", "dark")
  h <- plot_panel_add_effect_times(h, true_teff, group_by, col1, 1, 0.7)
  h <- plot_panel_add_effect_times(h, signal_teff, group_by, col2, 2, 0.7)

  # Add signal
  aes <- plot_panel_create_aes(signal, NULL)
  h <- h + ggplot2::geom_line(
    data = signal,
    mapping = aes,
    inherit.aes = FALSE,
    color = linecolors[1]
  )

  # Add data
  h <- h + ggplot2::geom_point()

  return(h)
}

#' Initialize the panel plot
#'
#' @inheritParams plot_panel
#' @return a \code{ggplot} object
plot_panel_create <- function(data, i_test, y_transform, nrow, ncol) {

  # Create dummy variable for type (train/test)
  show_test <- !is.null(i_test)
  if (show_test) {
    n <- dim(data)[1]
    a <- rep("train", n)
    a[i_test] <- "test"
    type_name <- "Type"
    data[[type_name]] <- as.factor(a)
  } else {
    type_name <- NULL
  }

  # Create aes
  aes <- plot_panel_create_aes(data, type_name)
  h <- ggplot2::ggplot(data, aes)

  # Create faceting
  group_by <- colnames(data)[1]
  f <- stats::as.formula(paste("~", group_by))
  h <- h + ggplot2::facet_wrap(f, nrow = nrow, ncol = ncol)
  if (show_test) {
    h <- h + ggplot2::scale_shape_manual(values = c(16, 16))
    col3 <- color_palette(3)[3]
    cols <- c(col3, "black")
    h <- h + ggplot2::scale_color_manual(values = cols)
  }
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
    h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "z"),
      na.rm = TRUE,
      data = df,
      color = linecolor,
      linetype = linetype,
      lwd = lwd
    )
  }
  return(h)
}
