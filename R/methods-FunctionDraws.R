#' @describeIn FunctionDraws Print a summary about the object.
setMethod("show", "FunctionDraws", function(object) {
  desc <- class_info("FunctionDraws")
  cat(desc)
})

#' @describeIn FunctionDraws  Visualization.
#' Optional arguments (\code{...}) are passed to
#' \code{\link{plot_FunctionDraws}}.
#' @param x a \linkS4class{FunctionDraws} object to visualize
#' @param y unused argument
setMethod(
  "plot",
  signature = c("FunctionDraws", "missing"),
  function(x, y, ...) {
    plot_FunctionDraws(x, ...)
  }
)

#' Not implemented
#'
#' @param object object to plot
plot_FunctionDraws <- function(object) {
  stop("not implemented")
}
