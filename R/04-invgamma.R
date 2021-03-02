
#' Density and quantile functions of the inverse gamma distribution
#'
#' @description Using the same parametrization as Stan. More info
#' \href{https://mc-stan.org/docs/2_24/functions-reference/inverse-gamma-distribution.html}{here}.
#' @param alpha positive real number
#' @param beta positive real number
#' @param x point where to compute the density
#' @param log is log-scale used?
#' @return density/quantile value
#' @name dinvgamma_stanlike
#' @family functions related to the inverse-gamma distribution

#' @rdname dinvgamma_stanlike
dinvgamma_stanlike <- function(x, alpha, beta, log = FALSE) {
  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  if (beta <= 0) {
    stop("beta must be positive")
  }
  t1 <- alpha * log(beta) - lgamma(alpha)
  t2 <- -1 * (alpha + 1) * log(x)
  t3 <- -beta / x
  log_p <- t1 + t2 + t3
  if (log) {
    return(log_p)
  } else {
    return(exp(log_p))
  }
}

#' @param p quantile (must be between 0 and 1)
#' @rdname dinvgamma_stanlike
qinvgamma_stanlike <- function(p, alpha, beta) {
  check_positive(alpha)
  check_positive(beta)
  check_interval(p, 0, 1)
  r <- stats::qgamma(1 - p, shape = alpha, rate = beta)
  return(1 / r)
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
