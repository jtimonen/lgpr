#' Finalize the lgpfit object after sampling
#'
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @param threshold Covariate selection threshold.
#' @param verbose Should some output be printed?
#' @return An updated object of class \code{lgpfit}.
postproc <- function(fit, threshold = 0.95, verbose = FALSE){
  
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
  FFF      <- get_function_component_samples(fit)
  n_smp    <- dim(FFF)[1]
  n_dat    <- dim(FFF)[2]
  n_cmp    <- dim(FFF)[3] - 2
  
  # Get average inferred components
  FFF_avg  <- apply(FFF, c(2,3), mean)
  FFF_avg  <- matrix_to_df(FFF_avg)
  if(n != n_dat){ stop("Data size sanity check failed!") }
  
  
  # Compute relevances for each sample
  if(verbose){ cat("* Computing relevances over", n_smp, "posterior samples.\n") }
  p_comp   <- matrix(0, n_smp, n_cmp + 1)
  p_signal <- rep(0, n_smp)
  svar     <- rep(0, n_smp)
  evar     <- rep(0, n_smp)
  for(i in 1:n_smp){
    FFF_i           <- data.frame(FFF[i,,])
    colnames(FFF_i) <- colnames(FFF_avg)
    REL             <- compute_relevances(FFF_i, y_data, info, D)
    p_comp[i,]      <- REL$p_comp
    p_signal[i]     <- REL$p_signal
    svar[i]         <- REL$svar
    evar[i]         <- REL$evar
  }
  colnames(p_comp)  <- NAMES1
  
  if(verbose){ cat("\n") }
  
  # Set slot values
  fit@components         <- FFF_avg
  fit@relevances         <- list(samples = p_comp, 
                                 average = colMeans(p_comp))
  fit@signal_variance    <- svar
  fit@residual_variance  <- evar
  fit@selection          <- selection(fit, threshold, verbose)
  return(fit)
  
}


#' Select relevant components
#'
#' @export
#' @param object An object of class \code{lgpfit}.
#' @param threshold A threshold for proportion of explained variance
#' @param verbose should this print some output
#' @return the selected covariates
selection <- function(object, threshold = 0.95, verbose = FALSE)
{
  if(class(object)!="lgpfit") stop("Class of 'object' must be 'lgpfit'!")
  if(threshold > 1 || threshold < 0) stop("'threshold' must be between 0 and 1!")
  info       <- object@model@info
  rel        <- object@relevances$average
  rel        <- sort(rel, decreasing = TRUE)
  i_noise    <- which(names(rel)=="noise")
  rel_rem    <- rel[-i_noise]
  rel        <- c(rel[i_noise], rel_rem)
  names(rel) <- c("noise", names(rel_rem))
  for(j in 1:length(rel)){
    h <- sum(rel[1:j])
    if(h > threshold){
      selected <- names(rel[1:j])
      res <- list(selected = selected, ev_sum = h, 
                  prop_ev = rel, threshold = threshold)
      return(res)
    }
  }
  res <- list(selected = names(rel), ev_sum = 1,
              prop_ev = rel, threshold = threshold)
  return(res)
}


#' Assess convergence of the chains
#'
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @param skip_F Should F_mean, F_var etc. be ignored
#' @return A data frame with columns \code{c("Rhat", "Bulk_ESS", "Tail_ESS")}.
assess_convergence <- function(fit, skip_F = TRUE){
  m  <- rstan::monitor(fit@stan_fit, print = FALSE)
  m  <- as.data.frame(m)
  m  <- m[c("Rhat", "Bulk_ESS", "Tail_ESS")]
  if(skip_F){
    i1 <- which(grepl("F_mean_", rownames(m)))
    i2 <- which(grepl("F_var_", rownames(m)))
    m  <- m[-c(i1,i2),]
  }
  return(m)
}

