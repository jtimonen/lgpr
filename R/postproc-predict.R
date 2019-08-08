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


#' Preprocess some things before computing predictions
#'
#' @description This is a helper function for \code{\link{lgp_predict}}.
#' @param fit An object of class \code{lgpfit}.
#' @param X_test The test points where the predictions should be computed.
#' @param samples The samples argument to \code{\link{lgp_predict}}
predict_preproc <- function(fit, X_test, samples){
  
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  if(class(model)!="lgpmodel") stop("Class of 'fit@model' must be 'lgpmodel'!")
  minfo    <- model@info
  cnames   <- minfo$component_names
  stan_dat <- model@stan_dat
  n        <- stan_dat$n
  X_sd     <- t(stan_dat$X)
  X_sd     <- X_sd[1:n,]
  Xnn_sd   <- t(stan_dat$X_notnan)
  Xnn_sd   <- Xnn_sd[1:n] 
  D        <- stan_dat$D
  LH       <- stan_dat$LH
  
  # Other model info
  fields   <- c("HMGNS", "UNCRT", "caseID_to_rows", "row_to_caseID", 
                "DELTA", "USE_VAR_MASK", "vm_params", "cat_interact_kernel")
  info     <- stan_dat[fields]
  
  # Input checking
  if(LH!=1){
    stop("Computing predictions outside the data is possible only for models",
         " with Gaussian likelihood.")
  }
  cn1  <- colnames(X_test)
  cn2  <- colnames(X_sd)
  if(length(cn1) != length(cn2)){
    cn_str <- paste(colnames(X_sd), collapse=", ")
    stop("X_test must be a data frame with column names {", cn_str,
         "} (in this order).", sep="")
  }
  if(!all(cn1==cn2)){
    cn_str <- paste(colnames(X_sd), collapse=", ")
    stop("X_test must be a data frame with column names {", cn_str, 
         "} (in this order).", sep="")
  }
  
  # Edit X_test so that we are working the same scale as with X_sd
  # (continuous covariates scaled to have zero mean and unit variance)
  TSCL <- model@scalings$TSCL
  X_test[,2] <- TSCL$fun(X_test[,2])
  if(D[3]==1){
    inan <- which(Xnn_sd==0)
    X_sd[inan,3] <- NaN
  }
  
  # Scale also other continuous variables to same scale
  if(D[4]>0){
    CSCL <- model@scalings$CSCL
    for(j in 1:D[4]){
      cscl <- CSCL[[j]]
      ix <- 2 + D[3] + j
      X_test[,ix] <- cscl$fun(X_test[,ix])
    }
  }
  
  # X_data and y_data
  X_data <- X_sd
  X_data <- cbind(X_data, stan_dat$row_to_caseID)
  colnames(X_data)[dim(X_data)[2]] <- "caseID"
  
  X_test <- add_test_caseIDs(X_test, X_data)
  y_data <- as.numeric(stan_dat$y)
  
  # All option to sample indicides
  if(is.character(samples)){
    if(samples=="all"){
      PAR <- hyperparam_samples(fit)
      samples <- 1:dim(PAR)[1]
    }
  }
  
  info$case_ids <- get_case_ids(fit)
  
  # Return
  ret <- list(X_data   = X_data, 
              X_test   = X_test,
              y_data   = y_data, 
              samples  = samples,
              D        = D,
              info     = info,
              cnames   = cnames)
  return(ret)
}


#' Average predictions over samples
#'
#' @param LIST a list over samples
#' @return a list
average_predictions <- function(LIST){
  L      <- length(LIST)
  p1     <- LIST[[1]]
  mu_cmp <- 1/L * p1$mu_cmp
  s2_cmp <- 1/L * p1$s2_cmp
  mu_f   <- 1/L * p1$mu_f
  s2_f   <- 1/L * p1$s2_f
  s2_y   <- 1/L * p1$s2_y
  for(i in 2:L){
    p      <- LIST[[i]]
    mu_cmp <- mu_cmp + 1/L * p$mu_cmp
    s2_cmp <- s2_cmp + 1/L * p$s2_cmp
    mu_f   <- mu_f   + 1/L * p$mu_f
    s2_f   <- s2_f   + 1/L * p$s2_f
    s2_y   <- s2_y   + 1/L * p$s2_y
  }
  ret <- list(mu_cmp = mu_cmp,
              s2_cmp = s2_cmp,
              mu_f   = mu_f,
              s2_f   = s2_f,
              s2_y   = s2_y)
  return(ret)
}


#' Add case IDs to test data frame
#'
#' @param X_test test data frame
#' @param X_data data frame
#' @return Updated X_test data frame.
add_test_caseIDs <- function(X_test, X_data){
  d        <- dim(X_data)[2]
  id_data  <- X_data[,c(1,d)]
  id_test  <- X_test[,1]
  uid      <- unique(id_data)
  L        <- length(id_test)
  cid_test <- rep(0, L)
  for(i in 1:L){
    idx <- which(uid[,1]==id_test[i])
    cid_test[i] <- uid[idx,2]
  }
  X_test <- cbind(X_test, cid_test)
  colnames(X_test)[d] <- "caseID" 
  return(X_test)
}

#' Compute component-wise predictions at test points
#'
#' @description Used by \code{\link{lgp_predict}}.
#' @param X_data Covariate matrix (data points).
#' @param y_data Response variable (data points).
#' @param X_test Covariate matrix (test points).
#' @param params Kernel function and other hyperparameters
#' @param info other model info
#' @param D a vector of length 6
#' @param cnames Names of the model components.
#' @param TSCL time scaling function and its inverse
#' @param handle_extra What to do if test data contains individuals that are
#' not in the training data? Must be 'silent', 'warning' or 'error'.
#' @return A list.
compute_predictions <- function(X_data, y_data, X_test, 
                                params, D, info, cnames,
                                TSCL,
                                handle_extra = "warning"){
  
  # Check that the IDs are same
  uid1 <- unique(X_data[,1])
  uid2 <- unique(X_test[,1])
  extra <- setdiff(uid2, uid1)
  if(length(extra)>0){
    msg <- "The test data contains individuals that are not in the training data! {"
    msg <- paste(msg, paste(extra, collapse = ", "), "}", sep = "")
    if(handle_extra=="warning"){
      warning(msg)
    }else if(handle_extra=="error"){
      stop(msg)
    }else if(handle_extra=="silent"){
      #do nothing
    }else{
      stop("unknown keyword for handle_extra (", handle_extra, ")")
    }
  }
  
  # Get kernel hyperparameters and other needed parameters
  nam   <- names(params)
  alpha <- params[which(grepl("alpha_", nam))]
  ell   <- params[which(grepl("lengthscale_", nam))]
  if(D[3]){
    stp <- params[which(grepl("warp_steepness", nam))]
    if(info$HMGNS==0){
      beta <- params[which(grepl("beta", nam))]
    }else{
      beta <- NULL
    }
    if(info$UNCRT==1){
      t_ons <- params[which(grepl("T_onset", nam))]
      names(t_ons) <- info$case_ids
    }else{
      t_ons <- NULL
    }
  }else{
    stp  <- NULL
    beta <- NULL
    t_ons <- NULL
  }
  
  sigma_n <- params[which(grepl("sigma_n", nam))]
  KERNEL_INFO <- list(D       = D, 
                      alpha   = alpha, 
                      ell     = ell, 
                      stp     = stp, 
                      beta    = beta,
                      t_ons   = t_ons,
                      info    = info,
                      TSCL    = TSCL)
  
  # Compute kernel matrices
  KK   <- compute_kernel_matrices(X_data, X_data, KERNEL_INFO)
  KKs  <- compute_kernel_matrices(X_test, X_data, KERNEL_INFO)
  KKss <- compute_kernel_matrices(X_test, X_test, KERNEL_INFO)
  
  # Compute predictions
  PRED <- compute_predicted_components(KK, KKs, KKss, y_data, 
                                          sigma_n, info$DELTA)
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


#' Evaluate kernel matrices for each component
#'
#' @description Used by \code{\link{compute_predictions}}.
#' @param X1 Covariate matrix of size \code{n1} x \code{sum(D)}.
#' @param X2 Covariate matrix of size \code{n2} x \code{sum(D)}.
#' @param kernel_info A list of parameters and other kernel info.
#' @return An array of size \code{n1} x \code{n2} x \code{sum(D)}.
compute_kernel_matrices <- function(X1, X2, kernel_info)
{
  
  # Get hyperparams and other kernel info
  D     <- kernel_info$D
  alpha <- as.numeric(kernel_info$alpha)
  ell   <- as.numeric(kernel_info$ell)
  stp   <- as.numeric(kernel_info$stp)
  beta  <- as.numeric(kernel_info$beta)
  t_ons <- kernel_info$t_ons
  info  <- kernel_info$info
  VM    <- info$USE_VAR_MASK
  vm_params <- info$vm_params
  TSCL  <- kernel_info$TSCL
  cat_interact_kernel <- info$cat_interact_kernel
  
  # Get dimensions and age and id covariates
  n1   <- dim(X1)[1]
  n2   <- dim(X2)[1]
  d1   <- dim(X1)[2]
  d2   <- dim(X2)[2]
  if(d1!=d2){
    stop("matrices X1 and X2 must have the same number of columns!")
  }
  d    <- sum(D)
  KK   <- array(0, c(n1,n2,d))
  id1  <- X1[,1]
  id2  <- X2[,1]
  age1 <- X1[,2]
  age2 <- X2[,2]
  
  # Compute kernel matrices
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
    r        <- r + 1
    mag      <- alpha[r]
    len      <- ell[r]
    disAge1  <- X1[,3]
    disAge2  <- X2[,3]
    xnn1     <- as.numeric(!is.nan(disAge1))
    xnn2     <- as.numeric(!is.nan(disAge2))
    
    # alpha mask
    if(VM==1){
      K_var_mask <- compute_K_var_mask(disAge1, disAge2, vm_params, stp)
    }else{
      K_var_mask <- matrix(1, n1, n2)
    }
    
    # beta mask
    if(info$HMGNS==0){
      caseID_1 <- X1[,d1] # case ids must be the last column
      caseID_2 <- X2[,d2]
      K_beta   <- compute_K_beta(beta, caseID_1, caseID_2)
    }else{
      K_beta <- matrix(1, n1, n2)
    }
    
    # modify diseaseAges
    if(info$UNCRT==1){
      formatter    <- function(x){formatC(x, width = 2, format = "d", flag = "0")}
      id1_str  <- formatter(id1)
      id2_str  <- formatter(id2)
      t_onset1 <- t_ons[id1_str]
      t_onset2 <- t_ons[id2_str]
      disAge1  <- TSCL$fun_inv(age1) - t_onset1
      disAge2  <- TSCL$fun_inv(age2) - t_onset2
      na1 <- which(is.na(disAge1))
      na2 <- which(is.na(disAge2))
      disAge1[na1] <- NaN
      disAge2[na2] <- NaN
    }
    
    # Evaluate the nonstationary kernel
    KK[,,r] <- K_var_mask * K_beta * kernel_bin(xnn1, xnn2) * 
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
      if(cat_interact_kernel == 1){
        KK[,,r] <- kernel_cat(x1, x2) * kernel_se(age1, age2, mag, len)
      }else{
        KK[,,r] <- kernel_bin(x1, x2) * kernel_se(age1, age2, mag, len)
      }
      
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


#' Compute log-posterior predictive density at test points
#'
#' @export
#' @param PRED predictions
#' @param y_test values of the response variable at the test points
#' @return a matrix with size \code{n_samples} x \code{n_data}
compute_lppd <- function(PRED, y_test){
  LIST <- PRED$LIST
  L    <- length(LIST)
  n    <- length(y_test)
  
  # Compute lppd for all samples and points separately
  LPPD <- matrix(0, L, n)
  for(j in 1:L){
    pred <- LIST[[j]]
    mu_y <- pred$mu_f
    s2_y <- pred$s2_y
    for(i in 1:n){
      LPPD[j,i] <- log_gaussian_density(y_test[i], mu_y[i], s2_y[i])
    }
  }
  return(LPPD)
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