#' Helper function for plot_warp
#'
#' @param par_summary summary of warping parameter
#' @param dis_age the x-axis values
#' @param color_scheme name of \code{bayesplot} color scheme
#' @return a \code{ggplot} object
plot_warp_helper <- function(par_summary,
                             dis_age,
                             color_scheme = "brightblue") {

  # Get colors and quantiles
  scheme <- bayesplot::color_scheme_get(color_scheme)
  color_line <- dollar(scheme, "dark")
  color_inner <- dollar(scheme, "light_highlight")
  color_outer <- dollar(scheme, "light")
  w_50 <- cpp_warp_input(dis_age, a = par_summary[6])
  w_75 <- cpp_warp_input(dis_age, a = par_summary[7])
  w_25 <- cpp_warp_input(dis_age, a = par_summary[5])
  w_025 <- cpp_warp_input(dis_age, a = par_summary[4])
  w_975 <- cpp_warp_input(dis_age, a = par_summary[8])
  DF <- data.frame(cbind(dis_age, w_50, w_75, w_25, w_025, w_975))

  # Create ggplot object
  aes_line <- ggplot2::aes_string(x = "dis_age", y = "w_50")
  aes_outer <- ggplot2::aes_string(ymin = "w_025", ymax = "w_975")
  aes_inner <- ggplot2::aes_string(ymin = "w_25", ymax = "w_75")
  h <- ggplot2::ggplot(DF, aes_line) +
    ggplot2::geom_ribbon(aes_outer, fill = color_outer) +
    ggplot2::geom_ribbon(aes_inner, fill = color_inner) +
    ggplot2::geom_line(color = color_line)

  # Add titles and labels
  h <- h + ggplot2::labs(x = "Input", y = "Warped input")
  subt <- paste("Median steepness =", round(par_summary[6], 3))
  h <- h + ggplot2::ggtitle("Input-warping function", subtitle = subt)
  return(h)
}
