#' Compute component-wise predictions at test points
#'
#' @description Used by \code{\link{lgp_predict}}.
#' @param X_data Covariate matrix (data points).
#' @param y_data Response variable (data points).
#' @param X_test Covariate matrix (test points).
#' @param params Kernel function and other hyperparameters.
#' @param D A vector of length 6.
#' @param info Other model info.
#' @param cnames Names of the model components.
#' @return A list containing predicted means and variances.
compute_predictions <- function(X_data, y_data, X_test, params, D, info, cnames){
  
  nam    <- names(params)
  # Get kernel hyperparameters and other needed parameters
  alpha  <- params[which(grepl("alpha_", nam))]
  ell    <- params[which(grepl("lengthscale_", nam))]
  if(D[3]){
    stp <- params[which(grepl("warp_steepness", nam))]
    if(info$HMGNS==0){
      beta <- params[which(grepl("beta", nam))]
    }else{
      beta <- NULL
    }
  }else{
    stp  <- NULL
    beta <- NULL
  }
  sigma_n <- params[which(grepl("sigma_n", nam))]
  
  # Compute kernel matrices
  KK   <- compute_kernel_matrices(X_data, X_data, D, alpha, ell, stp, beta, info)
  KKs  <- compute_kernel_matrices(X_test, X_data, D, alpha, ell, stp, beta, info)
  KKss <- compute_kernel_matrices(X_test, X_test, D, alpha, ell, stp, beta, info)
  
  # Compute predictions
  PRED  <- compute_predicted_components(KK, KKs, KKss, y_data, sigma_n, info$DELTA)
  F_mu  <- PRED$F_mu
  F_var <- PRED$F_var
  colnames(F_mu)  <- c(cnames, "f")
  colnames(F_var) <- c(cnames, "f")
  ret <- list(F_mu = F_mu, F_var = F_var)
  return(ret)
}


#' Compute component-wise predictions at test points
#'
#' @description Used by \code{\link{compute_predictions}}.
#' @param KK Kernel matrices data vs. data.
#' @param KKs Kernel matrices test vs. data.
#' @param KKss Kernel matrices test vs. test.
#' @param y_data Response variable.
#' @param sigma_n Noise standard deviation parameter.
#' @param DELTA Diagonal jitter that ensures pos. def. kernel.
#' @return A list containing predicted means and variances.
compute_predicted_components <- function(KK, KKs, KKss, y_data, sigma_n, DELTA){
  
  n <- length(y_data)
  p <- dim(KKs)[1]
  d <- dim(KK)[3]
  
  K_sum   <- apply(KK, c(1,2), sum)
  Ks_sum  <- apply(KKs, c(1,2), sum)
  Kss_sum <- apply(KKss, c(1,2), sum)
  Ky      <- K_sum + sigma_n^2*diag(n) + DELTA*diag(n)
  
  F_mu  <- matrix(0, p , d + 1)
  F_var <- matrix(0, p , d + 1)
  for(j in 1:d){
    Ks  <- KKs[,,j] 
    Kss <- KKss[,,j]
    F_mu[,j]  <- Ks %*% solve(Ky, y_data)
    F_var[,j] <- diag(Kss - Ks %*% solve(Ky, t(Ks))) 
  }
  F_mu[,d+1]  <- Ks_sum %*% solve(Ky, y_data)
  F_var[,d+1] <- diag(Kss_sum - Ks_sum %*% solve(Ky, t(Ks_sum)))
  
  ret   <- list(F_mu = F_mu, F_var = F_var)
  return(ret)
}

#' Create a matrix of test points
#'
#' @export
#' @param object An object of class \code{lgpmodel} or \code{lgpfit}
#' @param t_test Test time points (will be same for each individual).
#' @return A data frame.
create_test_points <- function(object, t_test){
  
  if(class(object)!="lgpfit"){
    if(class(object)=="lgpmodel"){
      model <- object
    }else{
      stop("Class of 'object' must be 'lgpfit' or 'lgpmodel'!")
    }
  }else{
    model <- object@model
  }
  sdat     <- model@stan_dat
  if(sdat$n_test != 0){
    warning("The model already contains test points!")
  }
  D        <- model@stan_dat$D
  if(D[4] > 0){
    stop("There are non-time-related continuous covariates in the data. There is no",
         " automatic way of generating the covariate matrix at test points based",
         " on only the test time points!")
  }
  SCL      <- model@scalings$TSCL
  X        <- t(sdat$X)
  Xnn      <- sdat$X_notnan
  X_test   <- create_X_star(X, D, t_test, SCL, Xnn)
  return(X_test)
}

#' Compute log-posterior predictive density at test points
#'
#' @export
#' @param PRED predictions
#' @param test_data test data
#' @return A data frame.
compute_lppd <- function(PRED, test_data){
  print(PRED)
  print(test_data)
  return("Loiko morko sisaan? Loiko morko sisaan??")
}
