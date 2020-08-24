#' Plot colors to use
#'
#' @export
#' @description Plot commands which specify the \code{col} argument should use
#' colors that are obtained using this function.
#' @param main Color name. Must be a valid \code{bayesplot} color name.
#' @param variant Must be one of {"light", "light_highlight", "mid",
#' "mid_highlight", "dark", "dark_highlight"}.
#' @return A hex value of the color.
colorset <- function(main, variant = "mid") {
  scheme <- bayesplot::color_scheme_get(scheme = main)
  col <- scheme[[variant]]
  if (is.null(col)) {
    stop("Invalid color!")
  }
  return(col)
}
