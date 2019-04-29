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


#' Evaluate kernel matrices for each component
#'
#' @description Used by \code{\link{compute_predictions}}.
#' @param X1 Covariate matrix of size \code{n1} x \code{sum(D)}.
#' @param X2 Covariate matrix of size \code{n2} x \code{sum(D)}.
#' @param D A vector of length 6.
#' @param alpha Marginal standard deviation parameters.
#' @param ell Lengthscale parameters
#' @param stp Input warp steepness parameter.
#' @param beta Individual-specific disease effect magnitudes
#' @param info Other model info
#' @return An array of size \code{n1} x \code{n2} x \code{sum(D)}.
compute_kernel_matrices <- function(X1, X2, D, 
                                    alpha, ell, stp, beta, 
                                    info)
{
  n1   <- dim(X1)[1]
  n2   <- dim(X2)[1]
  d    <- sum(D)
  KK   <- array(0, c(n1,n2,d))
  id1  <- X1[,1]
  id2  <- X2[,1]
  age1 <- X1[,2]
  age2 <- X2[,2]
  r <- 0
  if(D[1]==1){
    r   <- r + 1
    mag <- alpha[1]
    len <- ell[1]
    KK[,,r] <- kernel_cat(id1, id2) * kernel_se(age1, age2, mag, len)
  }
  if(D[2]==1){
    r   <- r + 1
    mag <- alpha[r]
    len <- ell[r]
    KK[,,r] <- kernel_se(age1, age2, mag, len)
  }
  if(D[3]==1){
    # TODO :include the special feature UNCRT too
    r       <- r + 1
    mag     <- alpha[r]
    len     <- ell[r]
    disAge1 <- X1[,3]
    disAge2 <- X2[,3]
    xnn1    <- as.numeric(!is.nan(disAge1))
    xnn2    <- as.numeric(!is.nan(disAge2))
    row_to_caseID_1 <- get_case_row_mappings(xnn1, id1)$row_to_caseID
    row_to_caseID_2 <- get_case_row_mappings(xnn2, id2)$row_to_caseID
    if(info$HMGNS==0){
      K_beta <- compute_K_beta(beta, row_to_caseID_1, row_to_caseID_2)
    }else{
      K_beta <- matrix(1, n1, n2)
    }
    KK[,,r] <- K_beta * kernel_bin(xnn1, xnn2) * 
      kernel_ns(disAge1, disAge2, mag, len, a = stp, b = 0, c = 1)
  }
  if(D[4]>0){
    for(j in 1:D[4]){
      r   <- r + 1
      mag <- alpha[r]
      len <- ell[r]
      idx <- 2 + D[3] + j
      x1  <- X1[,idx]
      x2  <- X2[,idx]
      KK[,,r] <- kernel_se(x1, x2, mag, len)
    }
  }
  if(D[5]>0){
    for(j in 1:D[5]){
      r   <- r + 1
      mag <- alpha[r]
      len <- ell[r]
      idx <- 2 + D[3] + D[4] + j
      x1  <- X1[,idx]
      x2  <- X2[,idx]
      KK[,,r] <- kernel_cat(x1, x2) * kernel_se(age1, age2, mag, len)
    }
  }
  if(D[6]>0){
    for(j in 1:D[6]){
      r   <- r + 1
      mag <- alpha[r]
      idx <- 2 + D[3] + D[4] + D[5] + j
      x1  <- X1[,idx]
      x2  <- X2[,idx]
      KK[,,r] <- mag^2 * kernel_cat(x1, x2)
    }
  }
  return(KK)
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


#' Create X_star
#'
#' @param X covariate matrix
#' @param D covariate type information
#' @param t_test Test time points (will be same for each individual).
#' @param SCL time scaling function and its inverse
#' @param X_notnan indicates where \code{X_diseaseAge} is not \code{NaN}
#' @return A data frame.
create_X_star <- function(X, D, t_test, SCL, X_notnan){
  
  X_id       <- X[,1]
  uid        <- unique(X_id)
  N          <- length(uid)
  n          <- length(X_id)
  d          <- dim(X)[2]
  X_age      <- X[,2]
  age_raw    <- SCL$fun_inv(X_age)
  p          <- length(t_test)
  X_star     <- matrix(0, N*p, d)
  if(D[3]){
    X_star[,3] <- rep(NaN, N*p)
  }
  Xnn_star   <- rep(0, N*p)
  X_star[,1] <- rep(uid, each = p)
  
  j <- 0
  
  for(id in uid){
    j <- j + 1
    i1 <- which(X_id==id)[1]
    inds <- ((j-1)*p + 1):(j*p)
    X_star[inds,2] <- t_test
    if(D[3]==1){
      if(X_notnan[i1]==1){
        t_onset <- age_raw[i1] - X[i1,3] 
        X_star[inds, 3] <- t_test - t_onset
        Xnn_star[inds] <- 1
      }
    }
    if(D[5]>0){
      for(k in 1:D[5]){
        X_star[inds, 2 + D[3] + k] <- X[i1, 2 + D[3] + k] 
      }
    }
    if(D[6]>0){
      for(k in 1:D[6]){
        X_star[inds, 2 + D[3] + D[5] + k] <- X[i1, 2 + D[3] + D[5] + k] 
      }
    }
  }
  
  X_star <- data.frame(X_star)
  colnames(X_star) <- colnames(X)
  return(X_star)
}
