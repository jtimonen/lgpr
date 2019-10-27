#' Compute a squared exponential kernel matrix
#'
#' @param x1 vector of length n
#' @param x2 vector of length m
#' @param alpha marginal std (default = 1)
#' @param ell lengthscale (default = 1)
#' @return A kernel matrix of size n x m
kernel_se <- function(x1,x2,alpha=1,ell=1){
  if(ell<=0){
    stop("ell must be positive!")
  }
  if(alpha <0){
    stop("alpha cannot be negative")
  }
  n1 <- length(x1)
  n2 <- length(x2)
  X1 <- matrix(rep(x1,each=n2),n1,n2,byrow=T)
  X2 <- matrix(rep(x2,n1),n1,n2,byrow=T)
  K  <- alpha^2*exp(-0.5*(X1-X2)^2/ell^2)
  return(K)
}


#' Compute a categorical kernel matrix
#'
#' @param x1 (integer) vector of length n
#' @param x2 (integer) vector of length m
#' @param alpha marginal std (default = 1)
#' @return A (binary) kernel matrix of size n x m
kernel_cat <- function(x1,x2,alpha=1){
  if(alpha <0){
    stop("alpha cannot be negative")
  }
  n1 <- length(x1)
  n2 <- length(x2)
  X1 <- matrix(rep(x1,each=n2),n1,n2,byrow=T)
  X2 <- matrix(rep(x2,n1),n1,n2,byrow=T)
  K  <- matrix(as.numeric(X1==X2),n1,n2)
  return(alpha^2*K)
}

#' Compute a zeros-sum kernel matrix
#'
#' @param x1 (integer) vector of length n
#' @param x2 (integer) vector of length m
#' @param M number of categories
#' @param alpha marginal std (default = 1)
#' @return A (binary) kernel matrix of size n x m
kernel_zerosum <- function(x1,x2,M,alpha=1){
  if(alpha <0){
    stop("alpha cannot be negative")
  }
  n1 <- length(x1)
  n2 <- length(x2)
  K  <- matrix(0, n1, n2)
  for(i in 1:n1){
    for(j in 1:n2){
      if(x1[i]==x2[j]){
        K[i,j] = 1.0;
      }else{
        K[i,j] = -1.0/(M-1); 
      }
    }
  }
  return(alpha^2*K)
}

#' Compute a binary kernel matrix
#'
#' @param x1 (integer) vector of length n
#' @param x2 (integer) vector of length m
#' @param alpha marginal std (default = 1)
#' @param pos_class the positive class label
#' @return A kernel matrix of size n x m
kernel_bin <- function(x1,x2=NULL,alpha=1,pos_class=1){
  if(alpha <0){
    stop("alpha cannot be negative")
  }
  n1 <- length(x1)
  n2 <- length(x2)
  X1 <- matrix(rep(x1,each=n2),n1,n2,byrow=T)
  X2 <- matrix(rep(x2,n1),n1,n2,byrow=T)
  K1 <- matrix(as.numeric(X1==pos_class),n1,n2)
  K2 <- matrix(as.numeric(X2==pos_class),n1,n2)
  return(alpha^2*K1*K2)
}

#' Compute a nonstationary kernel matrix using input warping
#'
#' @param x1 vector of length n
#' @param x2 vector of length m
#' @param alpha marginal std (default = 1)
#' @param ell lengthscale in the warped space
#' @param a steepness of the warping function rise
#' @param b location of the effective time window
#' @param c maximum range
#' @param nan_replace the value to use for replacing NaN values
#' @return A kernel matrix of size n x m
kernel_ns <- function(x1,x2=NULL,alpha=1,ell,a,b,c, nan_replace = 0){
  if(a <= 0){
    stop("a must be positive")
  }
  if(c <= 0){
    stop("c must be positive")
  }
  x1[is.nan(x1)] <- nan_replace
  x2[is.nan(x2)] <- nan_replace
  w1 <- warp_input(x1,a,b,c)
  w2 <- warp_input(x2,a,b,c)
  K  <- kernel_se(w1,w2,alpha,ell)
  return(K)
}

#' Warp inputs
#'
#' @param t a vector
#' @param a steepness of the rise
#' @param b location of the effective time window
#' @param c maximum range
#' @return a vector of warped inputs \code{w(t)}
warp_input <- function(t,a,b,c){
  w <- 2*c*(-0.5 + 1/(1+exp(-a*(t-b))) )
  return(w)
}


#' Compute the multiplier matrix K_beta (to eneable heterogeneous disease effect)
#'
#' @param beta a row vector of length \code{N_cases}
#' @param row_to_caseID_1 mapping from row index to case ID
#' @param row_to_caseID_2 mapping from row index to case ID
#' @return a matrix
compute_K_beta <- function(beta, row_to_caseID_1, row_to_caseID_2){
  n1   <- length(row_to_caseID_1)
  n2   <- length(row_to_caseID_2)
  BETA <- matrix(0, n1, n2)
  for(i in 1:n1){
    i_case <- row_to_caseID_1[i]
    if(i_case > 0){
      b1 <- beta[i_case]
    }else{
      b1 <- 0
    }
    for(j in 1:n2){
      j_case    <- row_to_caseID_2[j]
      if(j_case > 0){
        b2 <- beta[j_case]
      }else{
        b2 <- 0
      }
      BETA[i,j] <- sqrt(b1*b2)
    } 
  }
  return(BETA)
}


#' Compute the variance mask kernel matrix
#'
#' @param disAge1 disease-related age covariate vector of length \code{n1}
#' @param disAge2 disease-related age covariate vector of length \code{n2}
#' @param vm_params vector of two mask function parameters
#' @param stp input warping steepness
#' @param nan_replace value to replace nans in disAge vectors
#' @return a matrix of size \code{n1} x \code{n2}
compute_K_var_mask <- function(disAge1, disAge2, vm_params, stp, nan_replace = 0){
  disAge1[is.nan(disAge1)] <- nan_replace
  disAge2[is.nan(disAge2)] <- nan_replace
  mask_fun <- function(x,a){1/(1+exp(-a*x))}
  a  <- stp * vm_params[2];
  h  <- vm_params[1];
  if(h >= 1 || h <=0){
    stop("vm_params[1] must be between 0 and 1!")
  }
  if(vm_params[2] <= 0){
    stop("vm_params[2] must be greater than 0!")
  }
  
  r  <- 1/a*log(h/(1-h));
  s1 <- mask_fun(disAge1 - r, a)
  s2 <- mask_fun(disAge2 - r, a)
  M  <- tcrossprod(s1, s2)
  return(M)
}
