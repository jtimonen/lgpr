#' @describeIn GaussianPrediction Print a summary about the object.
setMethod("show", "GaussianPrediction", function(object) {
  comp_names <- names(object@f_comp_mean)
  D <- dim(object@f_comp_mean[[1]])
  desc <- class_info_fp("GaussianPrediction", comp_names, D)
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


#' @describeIn Prediction Print a summary about the object.
setMethod("show", "Prediction", function(object) {
  desc <- class_info_fp("Prediction")
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
