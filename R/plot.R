#' Visualize S4 class objects using the plot generic
#'
#' @param x does nothing
#' @param y does nothing
#' @param color_scheme bayesplot color scheme
#' @return a ggplot object
#' @name plot

#' @rdname plot
#' @param fit an object of class \linkS4class{lgpfit}
setMethod(
  f = "plot",
  signature = "lgpfit",
  definition = function(fit, x = 1, y = 1, color_scheme = "red") {
    h <- plot_relevances(fit, color_scheme = color_scheme)
    return(h)
  }
)
