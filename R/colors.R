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
#' @param n an integer from 1 to 6
#' @return an array of \code{n} hex values
color_palette <- function(n) {
  c1 <- colorset("brightblue", "mid_highlight")
  c2 <- colorset("red", "mid_highlight")
  c3 <- colorset("orange", "mid")
  c4 <- colorset("green", "mid_highlight")
  c5 <- colorset("gray", "dark")
  palette <- c(c1, c2, c3, c4, c5)
  if (n <= 5) {
    out <- palette[1:n]
  } else if (n <= 6) {
    palette <- unlist(bayesplot::color_scheme_get(scheme = "brightblue"))
    out <- as.character(palette)
  } else {
    stop("number of colors can be at most 6")
  }
  return(out)
}

#' Visualize a color palette
#'
#' @inheritParams color_palette
#' @return a \code{ggplot} object
plot_color_palette <- function(n) {
  colors <- color_palette(n)
  x <- rep(c(0, 1), n)
  y <- rep(c(1:n), each = 2)
  col <- as.factor(rep(colors, each = 2))
  df <- data.frame(x, y, col)
  aes <- ggplot2::aes_string(x = x, y = y, color = col, group = col)
  h <- ggplot2::ggplot(df) +
    ggplot2::geom_line(aes, lwd = 1)
  h <- h + ggplot2::scale_color_manual(values = colors)
  blank <- ggplot2::element_blank()
  h <- h + ggplot2::theme(
    axis.text = blank,
    axis.title = blank,
    axis.ticks = blank
  )
  h <- h + ggplot2::ggtitle("Colors")
  return(h)
}

#' Line alpha decider function
#'
#' @description A function that takes as parameter the number of lines and
#' returns a line alpha value (between 0 and 1).
#' @param num_lines number of lines to draw
#' @return a value between 0 and 1
line_alpha_fun <- function(num_lines) {
  alpha <- exp(-0.013 * num_lines)
  max(alpha, 0.01)
}
