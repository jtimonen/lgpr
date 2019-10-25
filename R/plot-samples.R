
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


#' Visualize posterior uncertainty in the disease effect times
#' 
#' @export
#' @description Can only be used if the uncertainty of effect time was modeled.
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
  ptitle <- "Posterior distribution of the inferred disease effect time"
  sd <- fit@model@stan_dat
  if(sd$UNCRT==0){
    stop("The disease effect time was not modeled as uncertain!")
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
  bayesplot::color_scheme_set(scheme = color_scheme)
  p     <- bayesplot::mcmc_dens(df)
  beta  <- "beta" # avoid note from R CMD check
  p <- p + ggplot2::xlab(expression(beta))
  
  str <- paste("Affected individuals: ", paste(names(which(aff)), collapse = ", "), sep="")
  form <- fit@model@info$formula
  p <- p + ggplot2::labs(
    subtitle = str,
    title = ptitle
  )
  return(p)
}

