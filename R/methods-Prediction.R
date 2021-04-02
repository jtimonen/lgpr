#' @describeIn Prediction Print a summary about the object.
setMethod("show", "Prediction", function(object) {
  desc <- class_info("Prediction")
  cat(desc)
})

#' @describeIn Prediction  Visualization.
#' Optional arguments (\code{...}) are passed to
#' \code{\link{plot_Prediction}}.
#' @param x a \linkS4class{Prediction} object to visualize
#' @param y unused argument
setMethod(
  "plot",
  signature = c("Prediction", "missing"),
  function(x, y, ...) {
    plot_Prediction(x, ...)
  }
)

#' Not implemented
#'
#' @param object object to plot
plot_Prediction <- function(object) {
  stop("not implemented")
}
