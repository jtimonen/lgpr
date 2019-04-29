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
#' @return A kernel matrix of size n x m
kernel_ns <- function(x1,x2=NULL,alpha=1,ell,a,b,c){
  if(a <= 0){
    stop("a must be positive")
  }
  if(c <= 0){
    stop("c must be positive")
  }
  x1[is.nan(x1)] <- 0
  x2[is.nan(x2)] <- 0
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
  BETA <- matrix(1, n1, n2)
  for(i in 1:n1){
    i_case <- row_to_caseID_1[i]
    if(i_case > 0){
      for(j in 1:n2){
        j_case    <- row_to_caseID_2[j]
        if(j_case > 0){
          b1        <- beta[i_case]
          b2        <- beta[j_case]
          BETA[i,j] <- sqrt(b1*b2)
        }
      } 
    }
  }
  return(BETA)
}
