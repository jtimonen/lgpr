#' Inverse gamma distribution pdf with the same parametrization
#' as in Stan
#'
#' @description See
#' \href{https://mc-stan.org/docs/2_21/functions-reference/inverse-gamma-distribution.html}{here}.
#' @param alpha positive real number
#' @param beta positive real number
#' @param x point where to compute the density
#' @param log should logarithm of the density be returned?
#' @return density or log-density at \code{x}
dinvgamma_stanlike <- function(x, alpha, beta, log = FALSE) {
  if (alpha <= 0) { stop("alpha must be positive") }
  if (beta <= 0) { stop("beta must be positive") }
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

#' Quantiles of inverse gamma distribution pdf with the same parametrization
#' as in Stan
#'
#' @description See
#' \href{https://mc-stan.org/docs/2_21/functions-reference/inverse-gamma-distribution.html}{here}.
#' @param alpha positive real number
#' @param beta positive real number
#' @param p quantile (must be between 0 and 1)
#' @return a positive real number
qinvgamma_stanlike <- function(p, alpha, beta) {
  r <- stats::qgamma(1 - p, shape = alpha, rate = beta)
  return(1 / r)
}


#' Plot pdf of inverse gamma distribution
#'
#' @export
#' @description Uses the same parametrization as described
#' \href{https://mc-stan.org/docs/2_21/functions-reference/inverse-gamma-distribution.html}{here}.
#' @param alpha positive real number
#' @param beta positive real number
#' @param by grid size
#' @param log if prior pdf should be on log scale
#' @param IQR inter-quantile range width
#' @param return_quantiles should this return a list
#' @return a ggplot object or a list
plot_invgamma <- function(alpha, beta, by = 0.01,
                          log = FALSE, IQR = 0.95,
                          return_quantiles = FALSE) {
  
  # Compute inter-quantile range
  if (IQR <= 0 || IQR >= 1) {
    stop("IQR must be on interval (0,1)")
  }
  delta <- (1 - IQR) / 2
  q1 <- qinvgamma_stanlike(delta, alpha = alpha, beta = beta)
  q2 <- qinvgamma_stanlike(1 - delta, alpha = alpha, beta = beta)

  max <- 1.2 * q2
  t <- seq(by, max, by = by)
  ylab <- "Prob."
  if (log) {
    ylab <- "Log prob."
  }
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
  line_color <- color_set('red')
  fill_color <- color_set('red_muted')
  p1 <- p1 + ggplot2::geom_area(
    data = df_area,
    mapping = ggplot2::aes(x = t, y = y),
    fill = fill_color
  )
  p1 <- p1 + ggplot2::geom_vline(xintercept = q1, color = line_color)
  p1 <- p1 + ggplot2::geom_vline(xintercept = q2, color = line_color)
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
