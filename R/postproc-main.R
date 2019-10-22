#' Finalize the lgpfit object after sampling
#'
#' @description Creates the \code{lgpfit} slots
#' \enumerate{
#'   \item \code{components} - Inferred components.
#'   \item \code{components_corrected} - Covariate-effect corrected components.
#'   \item \code{component_relevances} - Inferred component relevances.
#'   \item \code{covariate_relevances} - Inferred covariate relevances.
#'   \item \code{signal_variance} - Signal variance.
#'   \item \code{residual_variance} - Residual variance.
#'   \item \code{covariate_selection} -  Covariate selection info
#' }
#' all of which are lists that contain the fields \code{samples} and \code{average}.
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @param threshold Covariate selection threshold.
#' @param average_before_variance Should the variances be computed using average components?
#' @param sample_idx If supplied, this just returns the inferred components for one sample.
#' @param ell_smooth Defines how to determine smoothing lengthscale for corrected shared age 
#' effect inference. Possible options are
#' \enumerate{
#'    \item \code{"ell_shared"} (default) - the sampled lengthscale of the shared age 
#'    component is used as \code{ell_smooth}
#'    \item \code{"none"} - no correction will be performed
#'    \item A numeric argument that directly defines \code{ell_smooth}
#' }
#' @param ell_smooth_multip a multiplier for ell_smooth
#' @param verbose Should some output be printed?
#' @return An updated object of class \code{lgpfit}.
postproc <- function(fit, 
                     threshold               = 0.95, 
                     ell_smooth              = "ell_shared",
                     ell_smooth_multip       = 1,
                     sample_idx              = NULL,
                     average_before_variance = FALSE,
                     verbose                 = FALSE){
  
  model  <- fit@model
  info   <- model@info
  sdat   <- model@stan_dat
  y_data <- as.numeric(sdat$y)
  n      <- length(y_data)
  D      <- sdat$D
  n_cmp  <- length(info$component_names)
  NAMES1 <- c(info$component_names, "noise")
  NAMES2 <- c(info$covariate_names, "noise")
  
  # Get function component samples
  FFF      <- get_function_component_samples(fit, only_at_datapoints = TRUE)
  n_smp    <- dim(FFF)[1]
  n_dat    <- dim(FFF)[2]
  n_cmp    <- dim(FFF)[3] - 2
  
  # Get average inferred components
  FFF_avg  <- apply(FFF, c(2,3), mean)
  FFF_avg  <- matrix_to_df(FFF_avg)
  
  if(n != n_dat){
    stop("Data size sanity check failed!")
  }
  
  # Get shared age covariate and its lengthscale samples
  if(D[2]>0){
    ELL   <- rstan::extract(fit@stan_fit, pars = "lengthscale_sharedAge")[[1]]
    x_age <- as.numeric(sdat$X[2,])
  }else{
    ELL   <- NULL
    x_age <- NULL
  }
  
  if(is.null(sample_idx)){
    
    # Covariate and component relevance computations
    if(average_before_variance){
      if(verbose){ cat("* Averaging before computing variations.\n") }
    
      # Compute relevances using average F's
      ell          <- get_ell_smooth(ell_smooth, ell_smooth_multip, mean(ELL))
      REL          <- compute_relevances(FFF_avg, y_data, info, D, ell, x_age)
      FFF_cor_avg  <- REL$FFF_cor
      p_comp       <- REL$p_comp
      p_cov        <- REL$p_cov
      p_signal     <- REL$p_signal
      svar         <- REL$svar
      evar         <- REL$evar
      
    }else{
      
      # Compute relevances for each sample
      if(verbose){ cat("* Computing relevances over", n_smp, "posterior samples.\n") }
      p_comp   <- matrix(0, n_smp, n_cmp + 1)
      p_cov    <- matrix(0, n_smp, n_cmp + 1)
      p_signal <- rep(0, n_smp)
      svar     <- rep(0, n_smp)
      evar     <- rep(0, n_smp)
      FFF_cor_avg <- matrix(0, n, n_cmp + 2)
      for(i in 1:n_smp){
        ell             <- get_ell_smooth(ell_smooth, ell_smooth_multip, ELL[i])
        FFF_i           <- data.frame(FFF[i,,])
        colnames(FFF_i) <- colnames(FFF_avg)
        REL             <- compute_relevances(FFF_i, y_data, info, D, ell, x_age)
        p_comp[i,]      <- REL$p_comp
        p_cov[i,]       <- REL$p_cov
        p_signal[i]     <- REL$p_signal
        svar[i]         <- REL$svar
        evar[i]         <- REL$evar
        FFF_i_cor       <- REL$FFF_cor
        FFF_cor_avg     <- FFF_cor_avg + 1/n_smp*FFF_i_cor
      }
      colnames(FFF_cor_avg) <- colnames(FFF_i_cor)
      colnames(p_comp)      <- NAMES1
      colnames(p_cov)       <- NAMES2
    }
    if(verbose){ cat("\n") }
    
    # Set slot values
    fit@components           <- FFF_avg
    fit@components_corrected <- FFF_cor_avg
    fit@component_relevances <- list(samples = p_comp, average = colMeans(p_comp))
    fit@covariate_relevances <- list(samples = p_cov, average = colMeans(p_cov))
    fit@signal_variance      <- svar
    fit@residual_variance    <- evar
    fit@covariate_selection  <- varsel(fit, threshold, verbose)
    fit@postproc_info        <- list(threshold               = threshold,
                                     ell_smooth              = ell_smooth,
                                     ell_smooth_multip       = ell_smooth_multip,
                                     average_before_variance = average_before_variance)
    return(fit)
  }else{
    # Just get one sample
    if(average_before_variance){stop("Should not average if taking one sample!")}
    FFF_i           <- data.frame(FFF[sample_idx,,])
    colnames(FFF_i) <- colnames(FFF_avg)
    ell             <- get_ell_smooth(ell_smooth, ell_smooth_multip, ELL[sample_idx])
    if(verbose){ cat("ell_smooth =", ell, "\n") }
    REL             <- compute_relevances(FFF_i, y_data, info, D, ell, x_age)
    FFF_i_cor       <- REL$FFF_cor
    out             <- list(components = FFF_i, corrected = FFF_i_cor)
    return(out)
  }
  
}


#' Covariate selection
#'
#' @export
#' @param object An object of class \code{lgpfit}.
#' @param threshold A threshold for proportion of explained variance
#' @param verbose should this print some output
#' @return the selected covariates
varsel <- function(object, threshold = 0.95, verbose = FALSE)
{
  if(class(object)!="lgpfit") stop("Class of 'object' must be 'lgpfit'!")
  if(threshold > 1 || threshold < 0) stop("'threshold' must be between 0 and 1!")
  info       <- object@model@info
  rel        <- object@covariate_relevances$average
  rel        <- sort(rel, decreasing = TRUE)
  i_noise    <- which(names(rel)=="noise")
  rel_rem    <- rel[-i_noise]
  rel        <- c(rel[i_noise], rel_rem)
  names(rel) <- c("noise", names(rel_rem))
  for(j in 1:length(rel)){
    h <- sum(rel[1:j])
    if(h > threshold){
      selected <- names(rel[1:j])
      if(verbose){
        cat("* The following covariates explain ", round(h*100,2), "% of variance: {", sep="")
        str <- paste(selected, collapse=", ")
        cat(str)
        cat("}.\n")
      }
      return(list(selected = selected, ev_sum = h, prop_ev = rel))
    }
  }
  return(list(selected = names(rel), ev_sum = 1, prop_ev = rel))
}


#' Assess convergence of the chains
#'
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @param verbose should convergence info be printed?
#' @param recompute Should the Rhat statistics be recomputed?
#' @return Potential scale reduction factors (R_hat).
assess_convergence <- function(fit, verbose = FALSE, recompute = F){
  
  if(length(fit@Rhat)==0 || recompute){
    info     <- fit@model@info
    stan_fit <- fit@stan_fit
    LH       <- info$likelihood
    smr      <- rstan::summary(stan_fit)
    Rhat     <- smr$summary[,"Rhat"]
    
    # Skip  derived quantities
    if(!info$sample_F){
      i1 <- which(grepl("F_mean_cmp", names(Rhat)))
      i2 <- which(grepl("F_var_cmp", names(Rhat)))
      i3 <- which(grepl("F_mean_tot", names(Rhat)))
      i4 <- which(grepl("F_var_tot", names(Rhat)))
      i_skip <- unique(c(i1, i2, i3, i4))
    }else{
      i_skip <- which(grepl("F", names(Rhat)))
    }
    Rhat   <- Rhat[-i_skip]
  }else{
    Rhat <- fit@Rhat
  }
  
  rmax <- max(Rhat)
  imax <- which(Rhat == rmax)
  imax <- imax[1]
  rmax <- round(rmax, 3)
  if(verbose) cat("* The largest R_hat value (ignoring generated quantities) is ", rmax, " (", names(Rhat)[imax], ").\n", sep="")
  return(Rhat)
}

