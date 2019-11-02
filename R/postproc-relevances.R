#' Compute component relevances and estimate amount of noise
#' (one MCMC sample)
#'
#' @export
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @param relevance_method Component relevance determination method. 
#' Must be either \code{"f_mean"} or \code{"alpha"}.
#' @param verbose Should some output be printed?
#' @return An updated object of class \code{lgpfit}.
postproc_relevances <- function(fit, 
                                relevance_method = 'f_mean', 
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
    rel[i, ] <- compute_relevances(pars, model, relevance_method)
  }
  
  colnames(rel)  <- c(cnam, "noise")
  if(verbose){ cat("\n") }
  rel <- as.data.frame(rel)
  rel_avg <- colMeans(rel)
  rel_avg <- rel_avg/sum(rel_avg)
  res <- list(samples = rel, average = rel_avg)
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
#' @return a matrix of size \code{1} x \code{n_components + 1}
compute_relevances <- function(pars, model, method = 'f_mean'){
  if(method=='f_mean'){
    ret <- compute_relevances_fmean(pars, model)
  }else if(method=='alpha'){
    ret <- compute_relevances_alpha(pars, model)
  }else{
    stop("Invalid relevance method name (", method, ")!")
  }
  return(ret)
}


#' The f_mean relevance determination method 
#' 
#' @inheritParams compute_relevances
#' @return a matrix of size \code{1} x \code{n_components + 1}
compute_relevances_fmean <- function(pars, model){
  
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
  
  # Component relevances
  p_signal <- svar/(svar + evar)
  if(n_cmp > 1){
    p_comp  <- apply(FFF[,1:n_cmp], 2, stats::var)
    p_comp  <- p_comp/sum(p_comp)*p_signal  
  }else{
    p_comp <- p_signal
  }
  nam           <- c(model@info$component_names, "noise")
  p_comp        <- c(p_comp, 1-p_signal)
  names(p_comp) <- nam
  p_comp        <- t(as.matrix(p_comp))
  return(p_comp)
}


#' The alpha relevance determination method 
#' 
#' @inheritParams compute_relevances
#' @return a matrix of size \code{1} x \code{n_components + 1}
compute_relevances_alpha <- function(pars, model){
  
  stop("alpha relevance method not yet implemented")
  
  # Create the returned list
  ret <- 0#list(p_comp   = p_comp,
          #    p_signal = p_signal,
          #    svar     = svar,
          #    evar     = evar)
  return(ret)
}
