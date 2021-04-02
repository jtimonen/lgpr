#' @describeIn Prediction Print a summary about the object.
setMethod("show", "Prediction", function(object) {
  desc <- class_info("Prediction")
  cat(desc)
})

#' @describeIn Prediction  Visualization.
#' Optional arguments (\code{...}) are passed to
#' \code{\link{plot_pred}}.
#' @param x a \linkS4class{Prediction} object to visualize
#' @param y unused argument
setMethod(
  "plot",
  signature = c("Prediction", "missing"),
  function(x, y, ...) {
    plot_pred(pred = x, ...)
  }
)
