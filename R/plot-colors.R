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

#' A color palette function
#'
#' @param value an integer
color_palette <- function(value) {
  c1 <- colorset("brightblue", "mid_highlight")
  c2 <- colorset("red", "mid_highlight")
  c3 <- colorset("orange", "mid")
  c4 <- colorset("gray", "dark_highlight")
  palette <- c(c1, c2, c3, c4)
  if (value <= 4) {
    out <- palette[1:value]
  } else if (value <= 6) {
    palette <- unlist(bayesplot::color_scheme_get(scheme = "brightblue"))
    out <- as.character(palette[(7 - value):6])
  } else {
    stop("number of colors can be at most 6")
  }
  return(out)
}
