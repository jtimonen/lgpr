
#' Plot a simulated longitudinal data set for each individual separately
#'
#' @export
#' @param simData a list returned by \code{\link{simulate_data}}
#' @param linecolor line color
#' @param nrow an argument for \code{ggplot2::facet_wrap}
#' @param ncol an argument for \code{ggplot2::facet_wrap}
#' @param i_test test point indices
#' @param color_point data point color
#' @param color_test test point color
#' @param y_transform function to transform the data y
#' @param signal_name name of signal
#' @seealso For plotting each component separately,
#'  see \code{\link{plot_components_simdata}}
#' @return a ggplot object
plot_simdata <- function(simData,
                         linecolor = "gray50",
                         nrow = NULL,
                         ncol = NULL,
                         i_test = NULL,
                         color_point = "black",
                         color_test = "steelblue2",
                         signal_name = "signal",
                         y_transform = function(x) {
                           x
                         }) {
  vlinecolors <- c("firebrick3", "firebrick4")
  vlinetypes <- c(1, 5)
  dat <- simData$data
  comp <- simData$components
  g <- comp$g
  y <- y_transform(dat$y)
  yval <- c(g, y)
  leg <- rep(c(signal_name, "y (data)"), each = length(g))
  id <- rep(dat$id, 2)
  age <- rep(dat$age, 2)
  N <- length(unique(dat$id))
  n <- length(dat$id)

  DF <- data.frame(id, age, yval, leg)
  DF$age <- as.numeric(age)
  DF$yval <- as.numeric(yval)
  is_test <- rep("y_train", n)
  if (!is.null(i_test)) {
    is_test[i_test] <- "y_test"
  }
  it <- rep(is_test, 2)
  DF$is_test <- as.factor(it)

  # Create ggplot object
  h <- ggplot2::ggplot(data = DF, ggplot2::aes_string(
    x = "age", y = "yval", group = "leg"))
  subt <- paste(n, " data points, ", N, " individuals", sep = "")

  # Faceting
  h <- h + ggplot2::facet_wrap(. ~ id, nrow = nrow, ncol = ncol)

  # Plot real and observed disease onsets
  ons1 <- simData$onsets
  ons2 <- simData$onsets_observed
  no1 <- sum(!is.nan(ons1))
  if (no1 > 0) {
    subt <- paste(subt, ". Vertical lines are the real effect time (solid) ",
      "and observed disease onset / initiation time (dashed).",
      sep = ""
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

  # Plot signal line and data points
  h <- h + ggplot2::geom_line(ggplot2::aes_string(linetype = "leg"),
    color = linecolor,
    lwd = 1,
    alpha = 1
  ) +
    ggplot2::scale_shape_manual(values = c(NA, 16)) +
    ggplot2::scale_linetype_manual(values = c(1, 0))

  if (is.null(i_test)) {
    h <- h + ggplot2::geom_point(ggplot2::aes_string(shape = "leg"),
      na.rm = TRUE, color = color_point
    )
  } else {
    h <- h + ggplot2::geom_point(ggplot2::aes_string(
      shape = "leg",
      color = "is_test"
    ),
    na.rm = TRUE
    )
  }

  h <- h + ggplot2::labs(x = "Age", y = "y")

  # Theme and titles
  h <- h + ggplot2::ggtitle("Simulated data", subtitle = subt)
  h <- h + ggplot2::theme_bw()
  h <- h + ggplot2::theme(legend.title = ggplot2::element_blank())
  h <- h + ggplot2::theme(legend.position = "top")

  # Point color and type
  if (!is.null(i_test)) {
    h <- h + ggplot2::scale_color_manual(values = c(color_test, color_point))
  }
  return(h)
}


#' Visualize the components of a simulated data set
#' @export
#' @param simData simulated data object (list)
#' @param time_is_xvar is the time variable the x-axis variable
#' in all subplots?
#' @param marker point marker
#' @param ... additional arguments for \code{\link{plot_components}}
#' @return an object returned by \code{ggpubr::ggarrange} or a list
plot_components_simdata <- function(simData,
                                    time_is_xvar = TRUE,
                                    marker = 16, ...) {
  data <- simData$data
  FFF <- simData$components
  model <- full_model(data)

  # Create plot
  nnn <- dim(FFF)[1]
  ddd <- dim(FFF)[2]
  FFF_1 <- FFF[, 1:(ddd - 3)]
  FFF <- array(0, dim = c(1, nnn, ddd - 3))
  FFF[1, , ] <- as.matrix(FFF_1)

  h <- plot_components(FFF, NULL, model, time_is_xvar,
    marker = marker, ...
  )
  return(h)
}

#' Create a plotting data frame for ggplot
#'
#' @description A helper function for \code{plot_simdata_by_component}.
#' @param simData An object created using \code{simulate_data}.
#' @return a data frame
create_simdata_plot_df <- function(simData) {
  data <- simData$data
  comp <- simData$components
  id <- data$id
  age <- data$age
  n <- length(id)
  d <- dim(comp)[2] - 4
  CMP <- comp[, 1:(d + 1)]
  cn <- colnames(CMP)
  cn <- cn[1:(length(cn) - 1)]
  cn <- simdata_colnames_pretty(cn)
  cn <- c(cn, "f")
  id <- rep(id, d + 1)
  age <- rep(age, d + 1)
  component <- rep(cn, each = n)
  value <- as.numeric(as.matrix(CMP))
  DF <- data.frame(
    id = as.factor(id),
    age = as.numeric(age),
    component = as.factor(component),
    value = as.numeric(value)
  )
  return(DF)
}

#' Simulated data column names in a prettier form
#'
#' @param cn column names
#' @return names of model components
simdata_colnames_pretty <- function(cn) {
  cn <- gsub("\\.", ", ", cn)
  cn_seq <- seq_len(length(cn))
  cn <- paste("f_", cn_seq, "(", cn, ")", sep = "")
  return(cn)
}
