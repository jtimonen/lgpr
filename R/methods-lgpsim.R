#' Visualize an lgpsim object
#'
#' @export
#' @param simdata an object of class \linkS4class{lgpsim}
#' @param ... other arguments to \code{\link{sim_plot}}
#' @return a \code{ggplot} object
plot_simdata <- function(simdata, ...) {
  sim_plot(simdata, ...)
}

#' @rdname plot_simdata
#' @param x an object of class \linkS4class{lgpsim}
#' @param y not used
#' @return a \code{ggplot} object
setMethod(
  f = "plot",
  signature = signature(x = "lgpsim", y = "missing"),
  definition = function(x, ...) {
    plot_simdata(x, ...)
  }
)
