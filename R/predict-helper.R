

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
    row_to_caseID_1 <- get_case_row_mappings(xnn1, id1, only_R2C = TRUE)
    row_to_caseID_2 <- get_case_row_mappings(xnn2, id2, only_R2C = TRUE)
    #print(row_to_caseID_1)
    #print(row_to_caseID_2)
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
