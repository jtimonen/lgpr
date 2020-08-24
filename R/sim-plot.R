#' Plot a simulated longitudinal data set for each individual separately
#'
#'
#' @export
#' @inheritParams sim_plot_create
#' @inheritParams sim_plot_add_signal_and_data
#' @param nrow number of rows, an argument for \code{ggplot2::facet_wrap}
#' @param ncol number of columns, an argument for \code{ggplot2::facet_wrap}
#' @seealso For plotting each component separately,
#'  see \code{\link{sim_plot_components}}
#' @return a ggplot object
sim_plot <- function(simdata,
                     nrow = NULL,
                     ncol = NULL,
                     i_test = NULL,
                     linecolor = colorset("gray"),
                     dotcolor = "black",
                     testcolor = colorset("brightblue", "dark"),
                     signal_name = "signal",
                     y_transform = function(x) {
                       x
                     }) {

  # Create and add facet
  created <- sim_plot_create(simdata, i_test, y_transform, signal_name)
  h <- created$plot
  N <- created$N
  h <- h + ggplot2::facet_wrap(. ~ id, nrow = nrow, ncol = ncol)
  subtitle <- " "

  # Add possible effect times (vertical lines)
  et <- simdata@effect_times
  updated <- sim_plot_add_onsets(h, et, N)
  h <- updated$plot
  subtitle <- paste0(subtitle, updated$info)

  # Plot signal line and data points
  h <- sim_plot_add_signal_and_data(h, i_test, linecolor, dotcolor, testcolor)

  # Add axes, theme and titles
  h <- h + ggplot2::labs(x = "Age", y = "y")
  h <- h + ggplot2::ggtitle("Simulated data", subtitle = subtitle)
  h <- h + ggplot2::theme_bw()
  h <- h + ggplot2::theme(legend.title = ggplot2::element_blank())
  h <- h + ggplot2::theme(legend.position = "top")

  return(h)
}


#' Initialize a simulated data plot
#'
#' @description A helper function for \code{sim_plot}.
#' @param simdata an object of class \linkS4class{lgpsim}
#' @param i_test indices of possible test points
#' @param y_transform function to transform the data y
#' @param signal_name name of signal in the legend
#' @return updated \code{ggplot} object
sim_plot_create <- function(simdata, i_test, y_transform, signal_name) {
  dat <- simdata@data
  comp <- simdata@components
  g <- comp$g
  y <- y_transform(dat$y)
  yval <- c(g, y)
  leg <- rep(c(signal_name, "y (data)"), each = length(g))
  id <- rep(dat$id, 2)
  age <- rep(dat$age, 2)
  n <- length(dat$id)
  N <- length(unique(dat$id))

  DF <- data.frame(id, age, yval, leg)
  DF$age <- as.numeric(age)
  DF$yval <- as.numeric(yval)
  is_test <- rep("y_train", n)
  if (!is.null(i_test)) {
    is_test[i_test] <- "y_test"
  }
  it <- rep(is_test, 2)
  DF$is_test <- as.factor(it)

  # Create and return ggplot object
  aes <- ggplot2::aes_string(x = "age", y = "yval", group = "leg")
  h <- ggplot2::ggplot(data = DF, aes)
  list(plot = h, subtitle = "", N = N, n = n)
}


#' Add real and observed disease onsets to plot
#'
#' @description A helper function for \code{sim_plot}.
#' @param h a \code{ggplot} object
#' @param effect_times a list
#' @param N number of individuals
#' @return a list containing the updated \code{ggplot} object and info
sim_plot_add_onsets <- function(h, effect_times, N) {
  vlinecolors <- c(colorset("red"), colorset("red", "dark"))
  vlinetypes <- c(1, 5)
  ons1 <- effect_times$true
  ons2 <- effect_times$observed
  no1 <- sum(!is.nan(ons1))
  info <- ""
  if (no1 > 0) {
    info <- paste0(
      "Vertical lines are the real effect time (solid) \n",
      "and observed disease onset / initiation time (dashed)."
    )
    vline1.data <- data.frame(z = ons1, id = c(1:N))
    vline2.data <- data.frame(z = ons2, id = c(1:N))

    # Plot real effect times
    h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "z"),
      na.rm = TRUE,
      data = vline1.data,
      color = vlinecolors[1],
      linetype = vlinetypes[1],
      lwd = 0.7
    )

    # Plot observed onsets
    h <- h + ggplot2::geom_vline(ggplot2::aes_string(xintercept = "z"),
      na.rm = TRUE,
      data = vline2.data,
      color = vlinecolors[2],
      linetype = vlinetypes[2],
      lwd = 0.7
    )
  }

  # Return
  list(plot = h, info = info)
}


#' Add generated signal and data to plot
#'
#' @description A helper function for \code{sim_plot}.
#' @param h a \code{ggplot} object
#' @inheritParams sim_plot_create
#' @param linecolor line color for signal
#' @param dotcolor dot color for data points
#' @param testcolor dot color for possible test points
#' @return the updated \code{ggplot} object
sim_plot_add_signal_and_data <- function(h, i_test,
                                         linecolor,
                                         dotcolor,
                                         testcolor) {

  # Plot signal line
  aes <- ggplot2::aes_string(linetype = "leg")
  h <- h + ggplot2::geom_line(aes, color = linecolor, lwd = 1, alpha = 1)
  h <- h + ggplot2::scale_shape_manual(values = c(NA, 16))
  h <- h + ggplot2::scale_linetype_manual(values = c(1, 0))

  # Plot data points
  if (is.null(i_test)) {
    aes <- ggplot2::aes_string(shape = "leg")
    h <- h + ggplot2::geom_point(aes, na.rm = TRUE, color = dotcolor)
  } else {
    aes <- ggplot2::aes_string(shape = "leg", color = "is_test")
    h <- h + ggplot2::geom_point(aes, na.rm = TRUE)
  }

  # Edit point types
  if (!is.null(i_test)) {
    h <- h + ggplot2::scale_color_manual(values = c(testcolor, dotcolor))
  }

  return(h)
}

#' Visualize the components of a simulated data set
#'
#' @export
#' @param simdata an object of class \linkS4class{lgpsim}
#' @param time_is_xvar is the time variable the x-axis variable
#' in all subplots?
#' @param marker point marker
#' @param ... additional arguments for \code{\link{plot_components}}
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
