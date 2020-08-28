#' Visualize a simulated longitudinal data set
#'
#' @param x an object of class \linkS4class{lgpsim}
#' @param y not used
#' @param ... keyword arguments passed to \code{\link{sim_plot}}
#' @return a \code{ggplot} object
setMethod(
  f = "plot",
  signature = signature(x = "lgpsim", y = "missing"),
  definition = function(x, ...) {
    sim_plot(x, ...)
  }
)

#' Visualize a simulated longitudinal data set
#'
#' @export
#' @description Creates plots where each observation unit has a separate panel.
#' @param simdata an object of class \linkS4class{lgpsim}
#' @param signal_name signal label to show in legend
#' @param signal_column name of the signal column in \code{simdata$components}
#' @param x_name name of x-axis variable
#' @param y_name name of y-axis variable
#' @param group_by grouping variable
#' @param ... keyword arguments to \code{\link{plot_panel}}
#' @return a \code{ggplot object}
sim_plot <- function(simdata, signal_name = "signal", signal_column = "f",
                     x_name = "age", y_name = "y", group_by = "id", ...) {
  data <- simdata@data
  df_points <- data[c(group_by, x_name, y_name)]
  df_lines <- df_points
  df_lines[[y_name]] <- simdata@components[[signal_column]]
  colnames(df_lines)[3] <- signal_name

  true_teff <- null_if_all_nan(simdata@effect_times$true)
  signal_teff <- null_if_all_nan(simdata@effect_times$observed)

  h <- plot_panel(
    data = df_points,
    signal = df_lines,
    true_teff = true_teff,
    signal_teff = signal_teff,
    ...
  )

  info <- paste0(
    "Vertical lines are the real effect time (solid) \n",
    "and observed disease onset / initiation time (dashed)."
  )
  h <- h + ggplot2::ggtitle("Simulated data", subtitle = info)
  return(h)
}

#' Visualize the components of a simulated data set
#'
#' @export
#' @param simdata an object of class \linkS4class{lgpsim}
#' @param time_is_xvar is the time variable the x-axis variable
#' in all subplots?
#' @param marker point marker
#' @param ... keyword arguments for \code{\link{plot_components}}
#' @return an object returned by \code{ggpubr::ggarrange} or a list
sim_plot_components <- function(simdata,
                                time_is_xvar = TRUE,
                                marker = 16, ...) {
  comp <- simdata@components
  nnn <- dim(comp)[1]
  ddd <- dim(comp)[2]
  comp_1 <- comp[, 1:(ddd - 3)]
  comp <- array(0, dim = c(1, nnn, ddd - 3))
  comp[1, , ] <- as.matrix(comp_1)
  plot_components(comp, NULL, time_is_xvar, marker = marker, ...)
}