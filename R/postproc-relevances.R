#' Compute component relevances and estimate amount of noise
#' (one MCMC sample)
#'
#' @export
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @param relevance_method Component relevance determination method. 
#' Must be either \code{"f_mean"} or \code{"alpha"}.
#' @param noise_method Noise level determination method. 
#' Currently must be \code{"SSE"}.
#' @param verbose Should some output be printed?
#' @return An updated object of class \code{lgpfit}.
postproc_relevances <- function(fit, 
                                relevance_method = 'f_mean',
                                noise_method = 'SSE',
                                verbose   = FALSE){
  model  <- fit@model
  cnam   <- model@info$component_names
  DF     <- as.data.frame(fit@stan_fit)
  n_smp  <- dim(DF)[1]
  n_cmp  <- length(cnam)
  
  # Compute relevances for each sample
  if(verbose){ 
    cat("* Computing relevances over", n_smp, "posterior samples.\n") 
  }
  rel <- matrix(0, n_smp, n_cmp + 1)
  for(i in 1:n_smp){
    pars     <- DF[i,]
    rel[i, ] <- compute_relevances(pars, 
                                   model, 
                                   relevance_method,
                                   noise_method)
  }
  
  colnames(rel)  <- c(cnam, "noise")
  if(verbose){ cat("\n") }
  rel <- as.data.frame(rel)
  rel_avg <- colMeans(rel)
  rel_avg <- rel_avg/sum(rel_avg)
  res <- list(samples = rel, 
              average = rel_avg,
              method = relevance_method,
              noise_method = noise_method)
  return(res)
}


#' Compute component relevances and estimate amount of noise
#' (one MCMC sample)
#' 
#' @param pars A data frame representing one parameter sample,
#' i.e one row of \code{as.data.frame(stanfit)}, where stanfit
#' is an object of class \code{stanfit}
#' @param model An object of class \code{lgpmodel}
#' @param method Relevance determination method. Must be either
#' \code{"f_mean"} or \code{"alpha"}.
#' @param noise_method Noise level estimation method.
#' @return a matrix of size \code{1} x \code{n_components + 1}
compute_relevances <- function(pars, 
                               model, 
                               method,
                               noise_method){
  
  # Estimate noise level
  p_noise  <- compute_noise_level(pars, model, noise_method)
  
  if(method=='f_mean'){
    p_comp <- compute_relevances_fmean(pars, model)
  }else if(method=='alpha'){
    p_comp <- compute_relevances_alpha(pars, model)
  }else{
    stop("Invalid relevance method name (", method, ")!")
  }
  TOL <- 1e-6
  err <- abs(1-sum(p_comp))
  if(err > TOL){stop('Sanity check 1 failed!')}
  nam           <- c(model@info$component_names, "noise")
  p_comp        <- (1 - p_noise) * p_comp
  p_comp        <- c(p_comp, p_noise)
  names(p_comp) <- nam
  p_comp        <- t(as.matrix(p_comp))
  err           <- abs(1-sum(p_comp))
  if(err > TOL){stop('Sanity check 2 failed!')}
  return(p_comp)
}


#' The f_mean relevance determination method 
#' 
#' @inheritParams compute_relevances
#' @return a vector of length \code{n_components}
compute_relevances_fmean <- function(pars, model){
  
  # Get components
  n_cmp    <- sum(model@stan_dat$D)
  FFF      <- get_function_components_from_df(pars, model)
  FFF      <- as.data.frame(FFF)
  
  # Component relevances
  if(n_cmp > 1){
    p_comp  <- apply(FFF[,1:n_cmp], 2, stats::var)
    p_comp  <- p_comp/sum(p_comp) 
  }else{
    p_comp <- 1
  }
  return(p_comp)
}


#' The alpha relevance determination method 
#' 
#' @inheritParams compute_relevances
#' @return a vector of length \code{n_components}
compute_relevances_alpha <- function(pars, model){
  n_cmp <- sum(model@stan_dat$D)
  if(n_cmp > 1){
    pn      <- names(pars)
    ia      <- which(grepl('alpha_', pn))
    alpha   <- pars[ia]
    p_comp  <- as.numeric(alpha^2)
    p_comp  <- p_comp/sum(p_comp) 
  }else{
    p_comp <- 1
  }
  return(p_comp)

}


#' Determine noise level
#' 
#' @inheritParams compute_relevances
#' @return a value between 0 and 1
compute_noise_level <- function(pars, model, noise_method){
  
  if(noise_method!='SSE'){
    stop('noise_method must be SSE, no other methods implemented!')
  }
  
  # Get components
  n_data   <- model@stan_dat$n
  n_cmp    <- sum(model@stan_dat$D)
  FFF      <- get_function_components_from_df(pars, model)
  FFF      <- as.data.frame(FFF)
  
  # Signal variance
  y_pred   <- FFF$g
  svar     <- stats::var(y_pred)
  
  # Noise variance
  y_data <- model@stan_dat$y
  resid  <- y_pred - y_data
  evar   <- sum(resid^2)/(n_data-1)
  p_signal <- svar/(svar + evar)
  return(1 - p_signal) 
}

