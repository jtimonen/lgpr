#' Package colors
#'
#' @export
#' @description Here is defined a unified set of colors to be used throughout
#' the package. Plot commands which specify the \code{col} argument should use a
#' color that is specified here.
#' @param name Color name.
#' @return A hex value of the color.
color_set <- function(name) {
  allowed <- c("blue", "red")
  colors <- c("#3399ff", "#b22222")
  idx <- argument_check(arg = name, allowed = allowed)
  return(colors[idx])
}
