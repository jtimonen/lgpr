#' Visualize input warping function with several steepness parameter values
#'
#' @param wrp a vector of values of the warping steepness parameter
#' @param x a vector of input values
#' @param color line color
#' @param alpha line alpha
#' @return a \code{ggplot} object
plot_inputwarp <- function(wrp,
                           x,
                           color = colorset("red", "dark"),
                           alpha = 0.5) {
  x <- sort(x)
  L <- length(x)
  S <- length(wrp)
  W <- matrix(0, S, L)
  for (i in seq_len(S)) {
    w <- cpp_warp_input(x, a = wrp[i])
    W[i, ] <- w
  }
  af <- as.factor(rep(1:S, each = L))
  df <- data.frame(rep(x, S), as.vector(t(W)), af)
  colnames(df) <- c("x", "w", "idx")

  # Create ggplot object
  aes <- ggplot2::aes_string(x = "x", y = "w", group = "idx")
  plt <- ggplot2::ggplot(df, aes)

  # Add titles and labels
  plt <- plt + ggplot2::labs(x = "Input", y = "Warped input") +
    ggplot2::ggtitle("Input-warping function")
  plt <- plt + ggplot2::ylim(-1.0, 1.0)

  # Plot the actual lines of interest
  plt <- plt + ggplot2::geom_line(color = color, alpha = alpha)
  return(plt)
}


#' Plot the inverse gamma-distribution pdf
#'
#' @inheritParams dinvgamma_stanlike
#' @param by grid size
#' @param IQR inter-quantile range width
#' @param return_quantiles should this return a list
#' @param linecolor line color
#' @param fillcolor fill color
#' @return a \code{ggplot} object
#' @family functions related to the inverse-gamma distribution
plot_invgamma <- function(alpha, beta, by = 0.01,
                          log = FALSE, IQR = 0.95,
                          return_quantiles = FALSE,
                          linecolor = colorset("red", "dark"),
                          fillcolor = colorset("red", "mid")) {

  # Compute inter-quantile range
  check_interval(IQR, 0, 1)
  delta <- (1 - IQR) / 2
  q1 <- qinvgamma_stanlike(delta, alpha = alpha, beta = beta)
  q2 <- qinvgamma_stanlike(1 - delta, alpha = alpha, beta = beta)

  max <- 1.2 * q2
  t <- seq(by, max, by = by)
  ylab <- if (log) "Log prob." else "Prob."
  y <- dinvgamma_stanlike(t, alpha = alpha, beta = beta, log = log)
  df <- data.frame(t, y)

  # Plot prior
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = y)) +
    ggplot2::xlab("Steepness") +
    ggplot2::ylab(ylab) +
    ggplot2::theme_bw()

  # Plot IQR
  t <- seq(q1, q2, length.out = length(t))
  y <- dinvgamma_stanlike(t, alpha = alpha, beta = beta, log = log)
  df_area <- data.frame(t, y)
  p1 <- p1 + ggplot2::geom_area(
    data = df_area,
    mapping = ggplot2::aes(x = t, y = y),
    fill = fillcolor
  )
  p1 <- p1 + ggplot2::geom_vline(xintercept = q1, color = linecolor)
  p1 <- p1 + ggplot2::geom_vline(xintercept = q2, color = linecolor)
  p1 <- p1 + ggplot2::geom_line()
  q1r <- round(q1, 3)
  q2r <- round(q2, 3)
  msg <- paste0("The ", IQR * 100, "% IQR is [", q1r, ", ", q2r, "].\n")
  main <- paste0("Inv-Gamma(", alpha, ", ", beta, ")")
  p1 <- p1 + ggplot2::ggtitle(main, subtitle = msg)
  if (return_quantiles) {
    out <- list(
      plot = p1, lower = q1, upper = q2,
      mean = beta / (alpha - 1), IQR = IQR
    )
  } else {
    out <- p1
  }
  return(out)
}
