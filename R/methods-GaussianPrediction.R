#' @describeIn GaussianPrediction Print a summary about the object.
setMethod("show", "GaussianPrediction", function(object) {
  desc <- class_info("GaussianPrediction")
  cat(desc)
})

#' @describeIn GaussianPrediction Visualization.
#' Optional arguments (\code{...}) are passed to
#' \code{\link{plot_pred}}.
#' @param x a \linkS4class{GaussianPrediction} object to visualize
#' @param y unused argument
setMethod(
  "plot",
  signature = c("GaussianPrediction", "missing"),
  function(x, y, ...) {
    plot_pred(pred = x, ...)
  }
)
