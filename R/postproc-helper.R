
#' Get values of sampled function components at data points
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @return An array of size \code{n_samples} x \code{n_data} x \code{n_components+2}
get_function_component_samples <- function(fit)
{
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  info  <- fit@model@info
  sf    <- fit@stan_fit
  sd    <- fit@model@stan_dat
  LH    <- sd$LH
  C_hat <- sd$C_hat
  
  # Get function components
  if(!info$sample_F){
    F_cmp <- rstan::extract(sf, pars = "F_mean_cmp")$F_mean_cmp[,1,,]
    n_dim <- length(dim(F_cmp))
    if(n_dim==3){
      # many components
      F_cmp <- aperm(F_cmp, c(1,3,2))
      F_tot <- rstan::extract(sf, pars = "F_mean_tot")$F_mean_tot[,1,]
    }else{
      # only one component
      F_tot <- F_cmp
    }
    
  }else{
    F_cmp <- rstan::extract(sf, pars = "F")$F[,1,,]
    n_dim <- length(dim(F_cmp))
    if(n_dim==3){
      # many components
      F_cmp <- aperm(F_cmp, c(1,3,2))
      F_tot <- apply(F_cmp, c(1,2), sum)
    }else{
      # only one component
      F_tot <- F_cmp
    }
    
  }
  
  # Get sampled signal
  if(LH==1 || LH==0){
    G_tot <- F_tot
  }else if(LH==2 || LH==3){
    G_tot <- exp(F_tot + C_hat)
  }else if(LH==4){
    G_tot <- sd$N_trials * 1/(1 + exp(-F_tot))
  }else{
    stop("Unknown likelihood!")
  }
  
  # Concatenate
  n_smp <- dim(F_cmp)[1]
  n_tot <- dim(F_cmp)[2]
  if(n_dim == 3){
    n_cmp <- dim(F_cmp)[3]
  }else{
    n_cmp <- 1
  }
  FFF   <- array(0, c(n_smp, n_tot, n_cmp + 2))
  FFF[,,1:n_cmp]  <- F_cmp
  FFF[,,n_cmp+1]  <- F_tot
  FFF[,,n_cmp+2]  <- G_tot
  
  nam <- info$component_names
  nam <- c(nam, "f", "g")
  dimnames(FFF)[[3]] <- nam
  return(FFF)
}


#' Covariate and component relevance calculations
#' @param FFF a data frame of size \code{n_data} x \code{n_components+2}
#' @param y_data (scaled) measurements of the response variable
#' @param info model info
#' @param D a vector of length 6
#' @return a list
compute_relevances <- function(FFF, y_data, info, D){
  
  # Signal variance
  y_pred   <- FFF$g
  svar     <- stats::var(y_pred)
  
  # Noise variance
  resid  <- y_pred - y_data
  n_data <- dim(FFF)[1]
  evar   <- sum(resid^2)/(n_data-1)
  
  # Component relevances
  n_cmp <- length(info$component_names)
  p_signal <- svar/(svar + evar)
  if(n_cmp > 1){
    p_comp  <- apply(FFF[,1:n_cmp], 2, stats::var)
    p_comp  <- p_comp/sum(p_comp)*p_signal  
  }else{
    p_comp <- p_signal
  }
  nam           <- c(info$component_names, "noise")
  p_comp        <- c(p_comp, 1-p_signal)
  names(p_comp) <- nam
  p_comp        <- t(as.matrix(p_comp))
  
  # Create the returned list
  ret <- list(p_comp   = p_comp,
              p_signal = p_signal,
              svar     = svar,
              evar     = evar)
  return(ret)
}


#' A helper function 
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @return a list
get_predicted <- function(fit)
{
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  info <- fit@model@info
  sf   <- fit@stan_fit
  sdat <- fit@model@stan_dat
  n    <- sdat$n
  if(!info$sample_F){
    F_mean <- rstan::extract(sf, pars = "F_mean_tot")$F_mean_tot[,1,]
    F_var  <- rstan::extract(sf, pars = "F_var_tot")$F_var_tot[,1,]
    pred   <- colMeans(F_mean)
    std    <- sqrt(colMeans(F_var))
    ret    <- list(pred = pred, std = std)
  }else{
    F_smp <- rstan::extract(sf, pars = "F")$F[,1,,]
    F_smp <- aperm(F_smp, c(1,3,2))
    F_smp <- apply(F_smp, c(1,2), sum)
    LH    <- fit@model@stan_dat$LH
    C_hat <- fit@model@stan_dat$C_hat
    if(LH == 1){
      G_smp <- F_smp
    }else if(LH == 2 || LH == 3){
      G_smp <- exp(F_smp + C_hat)
    }else if(LH == 4){
      N_trials <- sdat$N_trials
      G_smp <- N_trials * 1/(1 + exp(-F_smp))
    }else{
      stop("Unknown likelihood!")
    }
    ret <- list(F_smp = F_smp, G_smp = G_smp)
  }
  return(ret)
}


#' Separate the covariate effects from an interaction components of a categorical 
#' covariate and age
#'
#' @param f_post a matrix of size \code{n} x \code{sum(D)}
#' @param D a vector of length 6
#' @param ell kernel lengthscale
#' @param t vector of \code{n} time points corresponding to \code{f_post}
#' @param i_edit Indices of columns whose effect should be moved to shared age.
#' @return a corrected \code{f_post}
separate_effects <- function(f_post, t, D, ell, i_edit){
  if(is.null(ell)){
    return(f_post)
  }else{
    i_shared <- D[1] + 1
    if(any(i_edit<=0)){stop("invalid input 'i_edit'!")}
    for(j in i_edit){
      fj                <- f_post[,j]
      ks                <- kernel_smoothing(fj, t, t, ell)
      f_post[,i_shared] <- f_post[,i_shared] + ks
      f_post[,j]        <- f_post[,j] - ks
    }
    return(f_post)
  }
}


#' Estimate conditional mean time profile using gaussian kernel smoothing
#'
#' @param v a vector of length \code{n} to be smoothed
#' @param t vector of \code{n} time points corresponding to \code{y}
#' @param t_out the set of \code{p} time points where the smoothing should be evaluated
#' @param ell kernel lengthscale
#' @return a vector of length \code{p}
kernel_smoothing <- function(v, t, t_out, ell){
  p <- length(t_out)
  n <- length(t)
  K <- kernel_se(t_out, t, alpha = 1, ell = ell)
  for(idx in 1:p){
    K[idx,] <- K[idx,]/sum(K[idx,]) # normalize kernel rows to sum to 1
  }
  y_out <- K%*%v
  return(y_out)
}


#' Get a posterior estimate of model (hyper)parameters
#'
#' @param object An (incomplete) object of class \code{lgpfit}.
#' @param type Must be "mean", "median", or "map".
#' @return a data frame
hyperparam_estimate <- function(object, type = "mean"){
  DF    <- as.data.frame(object@stan_fit)
  nam   <- names(DF)
  i_eta <- which(grepl("ETA",nam))
  i_f   <- which(grepl("F",nam))
  i_lp  <- which(grepl("lp__",nam))
  if(type=="mean"){
    DF <- colMeans(DF)
  }else if(type=="median"){
    DF <- apply(DF, 2, stats::median)
  }else if(type=="map"){
    i_max <- which(DF$lp__==max(DF$lp__))
    DF <- colMeans(DF[i_max,]) # colmeans just to get correct data type
  }else{
    stop("Invalid type!")
  }
  DF    <- DF[-c(i_f,i_eta,i_lp)]
  return(DF)
}


#' Get a set of model (hyper)parameter samples
#'
#' @param object An (incomplete) object of class \code{lgpfit}.
#' @param samples Sample indices. If NULL, all samples are taken.
#' @return a data frame
hyperparam_samples <- function(object, samples = NULL){
  nam   <- names(object@stan_fit)
  i_eta <- which(grepl("ETA",nam))
  i_f   <- which(grepl("F",nam))
  i_lp  <- which(grepl("lp__",nam))
  nam   <- nam[-c(i_f,i_eta,i_lp)]
  ext   <- rstan::extract(object@stan_fit, pars = nam)
  n     <- length(ext[[1]])
  cn    <- names(ext)
  d     <- length(cn)
  OUT   <- matrix(0, n, d)
  for(j in 1:d){
    OUT[,j] <- ext[[j]]
  }
  OUT <- data.frame(OUT)
  colnames(OUT) <- cn
  if(is.null(samples)){
    samples <- 1:n
  }
  OUT <- OUT[samples,]
  return(OUT)
}


#' Get average runtime of a chain
#' @export
#' @param object An object of class \code{lgpfit}.
#' @return Average runtimes for warmup and sampling
get_runtime <- function(object){
  TIM <- rstan::get_elapsed_time(object@stan_fit)
  n_chains <- dim(TIM)[1]
  t1 <- round(mean(TIM[,1]), 2)
  t2 <- round(mean(TIM[,2]), 2)
  return(list(warmup = t1, sampling = t2))
}


#' Select the affected individuals
#' @export
#' @param object An object of class \code{lgpfit}.
#' @param medians.return Should the medians of beta parameters also be returned?
#' @param threshold A value that the median of beta has to exceed 
#' @return A binary vector indicating the individuals for which the disease effect is inferred
#' to exist.
affected <- function(object, medians.return = FALSE, threshold = 0.5){
  DF <- as.data.frame(object@stan_fit)
  i_beta <- grep("beta", names(DF))
  if(length(i_beta)==0){
    stop("The disease effect was not modeled heterogeneously!")
  }
  BET      <- DF[i_beta]
  b50      <- apply(BET, 2, stats::median)
  affected <- (b50 >= threshold)
  cid      <- get_case_ids(object)
  names(affected) <- cid
  if(medians.return){
    return(list(affected = affected, medians = b50))
  }else{
    return(affected)
  }
}
