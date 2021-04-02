#' @describeIn GaussianPrediction Print a summary about the object.
setMethod("show", "GaussianPrediction", function(object) {
  desc <- class_info("GaussianPrediction")
  cat(desc)
})

#' @describeIn GaussianPrediction  Visualization.
#' Optional arguments (\code{...}) are passed to
#' \code{\link{plot_GaussianPrediction}}.
#' @param x a \linkS4class{GaussianPrediction} object to visualize
#' @param y unused argument
setMethod(
  "plot",
  signature = c("GaussianPrediction", "missing"),
  function(x, y, ...) {
    plot_GaussianPrediction(x, ...)
  }
)

#' Not implemented
#'
#' @param object object to plot
plot_GaussianPrediction <- function(object) {
  stop("not implemented")
}
