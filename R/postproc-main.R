#' Finalize the lgpfit object after sampling
#'
#' @inheritParams postproc_relevances
#' @inheritParams selection
#' @return An updated object of class \code{lgpfit}.
postproc <- function(fit, 
                     threshold = 0.95,
                     relevance_method = 'f_mean', 
                     verbose   = FALSE){
  
  # Set slot values
  fit@relevances <- postproc_relevances(fit, 
                                        relevance_method, 
                                        "SSE",
                                        verbose)
  fit@selection  <- selection(fit, threshold)
  return(fit)
}


#' Assess convergence of the chains
#'
#' @param fit An (incomplete) object of class \code{lgpfit}.
#' @param skip_F_gen Should F_mean, F_var etc. be ignored
#' @return A data frame with columns \code{c("Rhat", "Bulk_ESS", "Tail_ESS")}.
assess_convergence <- function(fit, skip_F_gen = TRUE){
  m  <- rstan::monitor(fit@stan_fit, print = FALSE)
  m  <- as.data.frame(m)
  m  <- m[c("Rhat", "Bulk_ESS", "Tail_ESS")]
  if(skip_F_gen){
    i1 <- which(grepl("F_mean_", rownames(m)))
    i2 <- which(grepl("F_var_", rownames(m)))
    if(length(i1)*length(i2)>0){
      m  <- m[-c(i1,i2),] 
    }
  }
  return(m)
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
