
#' Barplot of covariate relevances
#'
#' @export
#' @param object an object of class \code{lgpfit}
#' @param color_scheme bayesplot color scheme name
#' @return a ggplot object
plot_relevances <- function(object, color_scheme = "red")
{
  # Colors
  bpc <- bayesplot::color_scheme_get(color_scheme)
  fill  <- bpc$mid
  color <- bpc$mid_highlight
  
  # Covariates
  info      <- object@model@info
  Covariate <- c(info$covariate_names, "noise")
  Relevance <- as.vector(object@relevances$average)
  df <- data.frame(cbind(Relevance, Covariate))
  df$Relevance <- as.numeric(as.vector(df$Relevance))
  h <- ggplot2::ggplot(df, ggplot2::aes(x=Covariate, y=Relevance)) +
    ggplot2::geom_bar(stat = "identity", color = color, fill = fill) +
    ggplot2::theme_minimal()
  h <- h + ggplot2::labs(
    y = "Relevance",
    subtitle = paste("Model:", info$formula),
    title = "Relevances"
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
                           p            = 300,
                           color_scheme = "red",
                           b            = 0, 
                           c            = 1)
{
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  D <- model@stan_dat$D
  
  # Colors
  scheme      <- bayesplot::color_scheme_get(color_scheme)
  color_line  <- scheme$dark
  color_inner <- scheme$light_highlight
  color_outer <- scheme$light
  
  # Plot
  if(D[3]==1){
    X        <- t(model@stan_dat$X)
    X_disAge <- X[,3]
    ran      <- range(X_disAge)
    R        <- ran[2]- ran[1]
    ttt      <- seq(ran[1]-0.1*R, ran[2]+0.1*R, length.out = p) 
    
    sf    <- fit@stan_fit
    tsmr  <- rstan::summary(sf, pars = c("warp_steepness"))$summary
    w_50  <- warp_input(ttt, a = tsmr[6], b = b, c = c)
    w_75  <- warp_input(ttt, a = tsmr[7], b = b, c = c)
    w_25  <- warp_input(ttt, a = tsmr[5], b = b, c = c)
    w_025 <- warp_input(ttt, a = tsmr[4], b = b, c = c)
    w_975 <- warp_input(ttt, a = tsmr[8], b = b, c = c)
    
    diseaseAge <- ttt
    DF <- data.frame(cbind(diseaseAge, w_50, w_75, w_25, w_025, w_975))
    
    # Create ggplot object
    h <- ggplot2::ggplot(DF, ggplot2::aes_string(x = 'ttt', y = 'w_50')) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = 'w_025', ymax = 'w_975'),
                           fill = color_outer) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = 'w_25', ymax = 'w_75'),
                           fill = color_inner) +
      ggplot2::geom_line(color = color_line)
    
    h <- h + ggplot2::labs(x = "Disease-related age", y = "Warped input")
    subt <- paste("Median steepness =", round(tsmr[6], 3))
    h <- h + ggplot2::ggtitle("Input-warping function", subtitle = subt)
    return(h)
    
  }else{
    stop("Cannot visualize the input warping if 'diseaseAge' is not a model component.")
  }
}
