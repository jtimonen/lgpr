
#' Visualize the distribution of the model parameter samples
#' @export
#' @description This is a wrapper for functions in the \code{bayesplot} package.
#' @param object An object of class \code{lgpfit}.
#' @param pars parameter names
#' @param regex_pars regex for parameter names
#' @param point_est the point estimate type
#' @param transformations the parameter transformations
#' @param prob_outer outer interval
#' @param prob inner interval
#' @param binwidth width of histogram bins if \code{type = "hist"}
#' @param color_scheme See different color schemes in the \code{bayesplot} package.
#' @param off_diag_args Additional argument list for the pairs plot.
#' @param type Visualization type. Must be either \code{"dens"}, \code{"areas"},
#' \code{"intervals"}(default) or \code{"hist"}.
#' @param facet_args additional facetting arguments
#' @return a ggplot object
plot_samples <- function(object,
                         pars         = character(),
                         regex_pars   = character(),
                         type         = "intervals",
                         prob         = 0.5,
                         prob_outer   = 0.9,
                         color_scheme = "red",
                         point_est    = "median",
                         binwidth     = NULL,
                         transformations = list(),
                         off_diag_args   = list(size = 1),
                         facet_args   = list()
)
{
  if(class(object)!="lgpfit") stop("Class of 'object' must be 'lgpfit'!")
  posterior <- as.array(object@stan_fit)
  bayesplot::color_scheme_set(color_scheme)
  
  if(length(pars)==0 && length(regex_pars)==0){
    nam   <- names(object@stan_fit)
    i1    <- grep("F", nam)
    i2    <- grep("ETA", nam)
    iskip <- c(i1, i2, grep("lp__", nam))
    pars  <- nam[-iskip]
  }
  
  if(type=="areas"){
    h <- bayesplot::mcmc_areas(posterior,
                               pars       = pars,
                               regex_pars = regex_pars,
                               prob_outer = prob_outer,
                               prob       = prob,
                               point_est  = point_est,
                               transformations = transformations)
  }else if(type=="dens"){
    h <- bayesplot::mcmc_dens(posterior,
                              pars            = pars,
                              regex_pars      = regex_pars,
                              transformations = transformations,
                              facet_args      = facet_args)
  }else if(type=="intervals"){
    h <- bayesplot::mcmc_intervals(posterior,
                                   pars       = pars,
                                   regex_pars = regex_pars,
                                   prob_outer = prob_outer,
                                   prob       = prob,
                                   point_est  = point_est,
                                   transformations = transformations)
  }else if(type=="hist"){
    h <- bayesplot::mcmc_hist(posterior,
                              pars       = pars,
                              regex_pars = regex_pars,
                              binwidth   = binwidth,
                              transformations = transformations,
                              facet_args = facet_args)
  }else if(type=="pairs"){
    h <- bayesplot::mcmc_pairs(posterior, 
                               pars            = pars,
                               regex_pars      = regex_pars,
                               transformations = transformations,
                               off_diag_args   = off_diag_args)
  }else{
    stop(paste("Invalid type '", type, "'", sep=""))
  }
  return(h)
}


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
  Relevance <- as.vector(object@covariate_relevances$average)
  df <- data.frame(cbind(Relevance, Covariate))
  df$Relevance <- as.numeric(as.vector(df$Relevance))
  h <- ggplot2::ggplot(df, ggplot2::aes(x=Covariate, y=Relevance)) +
    ggplot2::geom_bar(stat = "identity", color = color, fill = fill) +
    ggplot2::theme_minimal()
  h <- h + ggplot2::labs(
    y = "Relevance",
    subtitle = paste("Model:", info$formula),
    title = "Covariate relevances"
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


#' Visualize posterior uncertainty in the disease onset
#' 
#' @export
#' @description Can only be used if the uncertainty of onset was modeled.
#' @param fit An object of class \code{lgpfit}.
#' @param color_scheme Name of bayesplot color scheme.
#' @param prob Inner interval
#' @param prob_outer Outer interval
#' @param point_est Point estimate type
#' @return a ggplot object
plot_onset <- function(fit, 
                       color_scheme = "red",
                       prob         = 1,
                       prob_outer   = 1,
                       point_est    = "none")
{
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  ptitle <- "Posterior distribution of the inferred disease onset"
  sd <- fit@model@stan_dat
  if(sd$UNCRT==0){
    stop("The disease onset was not modeled as uncertain!")
  }
  p <- plot_samples(fit, 
                    regex_pars   = "T_onset", 
                    type         = "areas", 
                    point_est    = point_est, 
                    prob         = prob, 
                    prob_outer   = prob_outer,
                    color_scheme = color_scheme)
  
  form <- fit@model@info$formula
  p <- p + ggplot2::labs(
    subtitle = paste("Model:", form),
    title = ptitle
  )
  return(p)
}


#' Visualize posterior samples of individual-specific disease effect magnitude
#' parameters
#' 
#' @export
#' @description Can only be used if the disease effect was modeled heterogeneously.
#' @param fit An object of class \code{lgpfit}.
#' @param color_scheme Name of bayesplot color scheme.
#' @param threshold Threshold for median.
#' @return a ggplot object
plot_beta <- function(fit, 
                      color_scheme = "red",
                      threshold    = 0.5)
{
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  ptitle <- "Posterior distribution of individual-specific disease effect magnitudes"
  
  aff   <- affected(fit, threshold = threshold)
  df    <- as.data.frame(fit@stan_fit)
  ibeta <- grep("beta", names(df))
  df    <- df[, ibeta]
  colnames(df) <- paste("id = ", names(aff), sep="")
  p <- bayesplot::mcmc_dens(df)
  beta <- "beta" # avoid note from R CMD check
  p <- p + ggplot2::xlab(expression(beta))
  
  str <- paste("Affected individuals: ", paste(names(which(aff)), collapse = ", "), sep="")
  form <- fit@model@info$formula
  p <- p + ggplot2::labs(
    subtitle = str,
    title = ptitle
  )
  return(p)
}

