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
#' @return A list.
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
  d  <- dim(F_mu)[2]
  
  # Predictions into return format
  mu_cmp <- F_mu[,1:(d-1)]
  s2_cmp <- F_var[,1:(d-1)]
  mu_f   <- F_mu[,d]
  s2_f   <- F_var[,d]
  s2_y   <- s2_f + sigma_n^2
  
  # Return
  ret <- list(mu_cmp = mu_cmp,
              s2_cmp = s2_cmp,
              mu_f   = mu_f,  
              s2_f   = s2_f,
              s2_y   = s2_y,
              pars   = params)
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
#' @param y_test values of the response variable at the test points
#' @param reduction must be either "mean", "sum" or "none"
#' @return a vector
compute_lppd <- function(PRED, y_test, reduction = "mean"){
  D  <- dim(PRED$F_mu)
  mu <- PRED$F_mu[, D[2]]
  s2 <- PRED$F_var[, D[2]]
  n  <- length(y_test)
  n1 <- length(mu)
  n2 <- length(s2)
  if(n!=n1 || n!=n2){
    stop("mu, sd and y_test must all have the same length! (", n, ", ", n1, ", ", n2, ")")
  }
  
  # Compute lppd for all points separately
  LPPD <- rep(0, n)
  for(i in 1:n){
    LPPD[i] <- log_gaussian_density(y_test[i], mu[i], s2[i])
  }
  
  # Return options
  if(reduction=="mean"){
    out <- mean(LPPD)
  }else if(reduction=="sum"){
    out <- sum(LPPD)
  }else if(reduction=="none"){
    out <- LPPD
  }else{
    stop("invalid value for reduction (", reduction, ")")
  }
  
  return(out)
}

#' Compute log-density for gaussian distribution
#'
#' @param x point x
#' @param mu mean
#' @param s2 variance
#' @return a number
log_gaussian_density <- function(x, mu, s2){
  p <- log(2*pi*s2) + 1/s2*(x-mu)^2
  p <- -0.5*p
  return(p)
}
