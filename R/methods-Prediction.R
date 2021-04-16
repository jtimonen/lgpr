#' @describeIn GaussianPrediction Print a summary about the object.
setMethod("show", "GaussianPrediction", function(object) {
  D1 <- num_components(object)
  D2 <- num_paramsets(object)
  D3 <- num_evalpoints(object)
  D <- c(D1, D2, D3)
  comp_names <- component_names(object)
  desc <- class_info_fp("GaussianPrediction", comp_names, D)
  cat(desc)
})

#' @describeIn Prediction Print a summary about the object.
setMethod("show", "Prediction", function(object) {
  D1 <- num_components(object)
  D2 <- num_paramsets(object)
  D3 <- num_evalpoints(object)
  D <- c(D1, D2, D3)
  comp_names <- component_names(object)
  desc <- class_info_fp("Prediction", comp_names, D)
  if (object@extrapolated) {
    origin_str <- paste0(
      "\nFunction draws have been extrapolated to x",
      " from original MCMC function draws at different x.",
      " Possibly also a reduction has been applied."
    )
    desc <- paste0(desc, origin_str)
  }
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

#' @describeIn GaussianPrediction Get names of components.
setMethod("component_names", "GaussianPrediction", function(object) {
  names(object@f_comp_mean)
})

#' @describeIn Prediction Get names of components.
setMethod("component_names", "Prediction", function(object) {
  names(object@f_comp)
})

#' @describeIn GaussianPrediction Get number of components.
setMethod("num_components", "GaussianPrediction", function(object) {
  length(component_names(object))
})

#' @describeIn Prediction Get number of components.
setMethod("num_components", "Prediction", function(object) {
  length(component_names(object))
})

#' @describeIn GaussianPrediction Get number of parameter combinations
#' (different parameter vectors) using which predictions were computed.
setMethod("num_paramsets", "GaussianPrediction", function(object) {
  D <- dim(object@f_comp_mean[[1]])
  D[1]
})

#' @describeIn Prediction Get number of parameter combinations
#' (different parameter vectors) using which predictions were computed.
setMethod("num_paramsets", "Prediction", function(object) {
  D <- dim(object@f_comp[[1]])
  D[1]
})

#' @describeIn GaussianPrediction Get number of points where
#' predictions were computed.
setMethod("num_evalpoints", "GaussianPrediction", function(object) {
  D <- dim(object@f_comp_mean[[1]])
  D[2]
})

#' @describeIn Prediction Get number of points where
#' predictions were computed.
setMethod("num_evalpoints", "Prediction", function(object) {
  D <- dim(object@f_comp[[1]])
  D[2]
})
