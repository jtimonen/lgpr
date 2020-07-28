
#' Barplot of covariate relevances
#'
#' @export
#' @param object an object of class \code{lgpfit}
#' @param violin Should a violin plot be used instead of a boxplot
#' @param color_scheme bayesplot color scheme name
#' @param ... Additional arguments to \code{ggplot2::geom_boxplot} or
#' \code{ggplot2::geom_violin}.
#' @return a ggplot object
plot_relevances <- function(object, violin = FALSE, color_scheme = "red", ...) {
  # Colors
  bpc <- bayesplot::color_scheme_get(color_scheme)
  fill <- bpc$mid
  color <- bpc$mid_highlight

  # Covariates
  info <- object@model@info
  Covariate <- c(info$component_names, "noise")
  Relevance <- as.matrix(object@relevances$samples)
  n_smp <- dim(Relevance)[1]
  n_cmp <- dim(Relevance)[2]
  rel <- as.numeric(Relevance)
  cname <- as.factor(rep(Covariate, each = n_smp))

  df <- data.frame(cname, rel)
  colnames(df) <- c("Component", "Relevance")
  h <- ggplot2::ggplot(df, ggplot2::aes_string(x = "Component", y = "Relevance"))
  if (violin) {
    h <- h + ggplot2::geom_violin(color = color, fill = fill, ...)
  } else {
    h <- h + ggplot2::geom_boxplot(color = color, fill = fill, ...)
  }

  h <- h + ggplot2::theme_minimal()
  h <- h + ggplot2::labs(
    y = "Relevance",
    title = paste0("Component relevances"),
    subtitle = paste0("Distribution over ", n_smp, " posterior samples")
  )
  return(h)
}


#' Visualize the input warping function for different parameter samples
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param color_scheme Name of bayesplot color scheme.
#' @param p number of plot points
#' @param b location of the effective time window (default = 0)
#' @param c maximum range (default = 1)
#' @return a ggplot object
plot_inputwarp <- function(fit,
                           p = 300,
                           color_scheme = "red",
                           b = 0,
                           c = 1) {
  if (class(fit) != "lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  D <- model@stan_dat$D

  # Colors
  scheme <- bayesplot::color_scheme_get(color_scheme)
  color_line <- scheme$dark
  color_inner <- scheme$light_highlight
  color_outer <- scheme$light

  # Plot
  if (D[3] == 1) {
    X <- t(model@stan_dat$X)
    X_disAge <- X[, 3]
    ran <- range(X_disAge)
    R <- ran[2] - ran[1]
    ttt <- seq(ran[1] - 0.1 * R, ran[2] + 0.1 * R, length.out = p)

    sf <- fit@stan_fit
    tsmr <- rstan::summary(sf, pars = c("warp_steepness"))$summary
    w_50 <- warp_input(ttt, a = tsmr[6], b = b, c = c)
    w_75 <- warp_input(ttt, a = tsmr[7], b = b, c = c)
    w_25 <- warp_input(ttt, a = tsmr[5], b = b, c = c)
    w_025 <- warp_input(ttt, a = tsmr[4], b = b, c = c)
    w_975 <- warp_input(ttt, a = tsmr[8], b = b, c = c)

    diseaseAge <- ttt
    DF <- data.frame(cbind(diseaseAge, w_50, w_75, w_25, w_025, w_975))

    # Create ggplot object
    h <- ggplot2::ggplot(DF, ggplot2::aes_string(x = "ttt", y = "w_50")) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "w_025", ymax = "w_975"),
        fill = color_outer
      ) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "w_25", ymax = "w_75"),
        fill = color_inner
      ) +
      ggplot2::geom_line(color = color_line)

    h <- h + ggplot2::labs(x = "Disease-related age", y = "Warped input")
    subt <- paste("Median steepness =", round(tsmr[6], 3))
    h <- h + ggplot2::ggtitle("Input-warping function", subtitle = subt)
    return(h)
  } else {
    stop("Cannot visualize the input warping if 'diseaseAge' is not a model component.")
  }
}
