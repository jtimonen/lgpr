#' Visualize S4 class objects using the plot generic
#'
#' @description Depending on the signature of the first argument, one
#' of the following signature-specific plotting functions is called
#' \itemize{
#'     \item \code{\link{plot_relevances}} for objects with signature
#'       \linkS4class{lgpfit}
#'     \item \code{\link{sim_plot}} for objects with signature
#'       \linkS4class{lgpsim}
#' }
#' @param x unused
#' @param y unused
#' @param ... keyword arguments to be passed to the signature-specific plotting
#' function
#' @return All methods return a \code{ggplot} object.
#' @name plot

#' @rdname plot
setMethod(
  f = "plot",
  signature = signature(x = "lgpfit", y = "missing"),
  definition = function(x, ...) {
    plot_relevances(x, ...)
  }
)

#' @rdname plot
setMethod(
  f = "plot",
  signature = signature(x = "lgpsim", y = "missing"),
  definition = function(x, ...) {
    sim_plot(x, ...)
  }
)
