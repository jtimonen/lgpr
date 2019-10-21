#' Create input for Stan
#'
#' @description Parses the \code{formula} and \code{data} input to \code{\link{lgp_model}}. Also performs
#' many input checks.
#' @inheritParams lgp_model
#' @param standardize Should the response variable be standardized?
#' @param varInfo Variable type info.
#' @param verbose Can this print some info?
#' @return A list containing the data to be given to \code{rstan::sampling}, some info about
#' preprocessing and all the information about scaling the inputs and response, and updated
#' variable type info.
#'
create_stan_input <- function(formula,
                              data,
                              prior,
                              likelihood,
                              varInfo,
                              standardize,
                              uncertain_effect_time,
                              equal_effect,
                              C_hat,
                              DELTA,
                              sample_F,
                              t_test,
                              verbose,
                              variance_mask,
                              cat_interact_kernel_type,
                              N_trials)
{
  
  # Check sample_F input
  if(!sample_F && likelihood!="Gaussian"){
    stop("sample_F must be true if likelihood is not Gaussian")
  }
  if(!sample_F && !is.null(t_test)){
    stop("t_test must be NULL if F is not sampled!")
  }
  
  # Check the model formula
  check_formula(formula, data)
  
  # Get the (scaled) response variable
  RESP       <- get_response(data, varInfo, standardize, likelihood)
  response   <- RESP$response
  
  # Resolve covariate types
  TYP      <- check_data(data, varInfo, verbose)
  types    <- TYP$types
  varInfo  <- TYP$varInfo
  
  # Take only those covariates that are needed by Stan
  covariates <- create_covariates_stan(data, varInfo, types, formula, verbose)
  D          <- covariates$D
  X          <- covariates$X
  X_notnan   <- covariates$X_notnan
  
  # Standardize continuous input covariates
  SCL <- standardize_inputs(X, D)
  X   <- SCL$X
  
  # Get some variables defining the model dimensions
  stan_dims <- get_model_dims(X, D, likelihood)
  
  # Set C_hat (constant)
  nb_or_pois <- likelihood %in% c("NB","Poisson")
  if(is.null(C_hat)){
    if(nb_or_pois){
      C_hat <- log(mean(response))
    }else{
      C_hat <- 0
    }
  }else{
    if(!nb_or_pois){
      stop("Only give the C_hat argument if observation model is Poisson or NB!")
    }
  }
  
  # Concatenate possible test points
  X_notnan <- as.vector(X_notnan)
  if(is.null(t_test)){
    X_final        <- X
    X_notnan_final <- X_notnan
    n_test         <- 0
  }else{
    X_star     <- create_X_star(X, D, t_test, SCL$TSCL, X_notnan)
    X_star[,2] <- SCL$TSCL$fun(X_star[,2])
    X_final    <- rbind(X, X_star)
    n_test     <- dim(X_star)[1]
    if(stan_dims$D[3]>0){
      Xnn_star       <- as.numeric(!is.nan(X[,3]))
      X_notnan_final <- c(X_notnan, Xnn_star)
    }else{
      X_notnan_final <- c(X_notnan, rep(0, n_test))
    }
  }
  
  # Categorical or binary kernel?
  if(cat_interact_kernel_type == "categorical"){
    cat_interact_kernel <- 1
  }else if(cat_interact_kernel_type == "binary"){
    cat_interact_kernel <- 0
  }else{
    stop("Invalid option '", cat_interact_kernel_type, ' for cat_interact_kernel_type!')
  }
  
  # Check N_trials
  if(is.null(N_trials)){
      N_trials <- rep(1, length(response))
  }else{
    if(likelihood != "binomial"){
      stop("Only give the N_trials argument if likelihood is binomial!")
    }
    if(length(N_trials)==1){
      N_trials <- rep(N_trials, length(response))
    }
    if(length(N_trials)!= length(response)){
      stop("Invalid length of N_trials!")
    }
  }

  
  # Create the list that is the Stan input
  stan_dat   <- list(X         = t(X_final),
                     X_id      = X_final[,1],
                     X_notnan  = X_notnan_final,
                     y         = response,
                     y_int     = round(response),
                     n_test    = n_test,
                     DELTA     = DELTA,
                     C_hat     = C_hat,
                     F_is_sampled = as.numeric(sample_F),
                     USE_VAR_MASK = as.numeric(variance_mask),
                     cat_interact_kernel = cat_interact_kernel,
                     N_trials  = N_trials
  )
  
  # Get some variables related to diseased individuals
  stan_dis  <- get_diseased_info(D, 
                                 X_final, 
                                 X_notnan_final,
                                 uncertain_effect_time, 
                                 equal_effect,
                                 SCL$TSCL)
  
  # Parse the prior
  stan_prior <- prior_to_stan(D, 
                              prior, 
                              stan_dis$HMGNS, 
                              stan_dis$UNCRT, 
                              stan_dis$N_cases,
                              stan_dis$T_observed,
                              stan_dis$T_last)
  
  # Return
  ret <- list(
    stan_dat  = c(stan_dat, stan_dims, stan_dis, stan_prior),
    varInfo   = varInfo,
    scalings  = list(TSCL = SCL$TSCL, CSCL = SCL$CSCL, YSCL = RESP$SCL)
  )
  return(ret)
  
}


#' Create the covariate matrix that is given to stan
#'
#' @param data the data frame that was passed to \code{lgp}
#' @param varInfo original variable type info
#' @param verbose can this print some info?
#' @param formula the model formula
#' @param types the types returned by \code{\link{check_data}}
#' @return a list 
create_covariates_stan <- function(data, varInfo, types, formula, verbose){
  
  # Create the design matrix X
  XD <- stan_input_X_and_D(data, varInfo, types, formula, verbose)
  
  # Do this because as.matrix() sometimes creates character matrices
  X <- data.matrix(XD$X) 
  D <- XD$D
  
  # Get location of NaNs because Stan does not accept NaNs
  if(D[3] > 0){
    X_notnan <- 1 - is.nan(X[,3])
    X[which(X_notnan==0),3] <- 0
  }else{
    X_notnan <- rep(0, length(X[,1]))
  }
  X_notnan <- as.vector(X_notnan)
  
  # If there are still NaNs, throw an error
  n_nans <- sum(is.nan(X))
  if(n_nans > 0){
    stop("Only the diseaseAge column of the data can contain NaNs!")
  }
  
  # Return
  covariates <- list(
    X        = X,
    X_notnan = X_notnan,
    D        = D
  )
  
  return(covariates)
}

