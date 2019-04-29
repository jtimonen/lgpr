#' The main function of the `lgpr` package
#'
#' @export
#' @description This is a wrapper for both \code{\link{lgp_model}} and \code{\link{lgp_fit}}.
#' It first creates an \code{lgpmodel} object and then fits the model, finally returning an
#' \code{lgpfit} object.
#' @inheritParams lgp_model
#' @inheritParams lgp_fit
#' @return An object of class \code{lgpfit}.
lgp <- function(formula,
                data,
                likelihood       = "Gaussian",
                prior            = prior_default(),
                uncertain_diagnosis = FALSE,
                equal_effect     = TRUE,
                id_variable      = "id",
                time_variable    = "age",
                disAge_variable  = NULL,
                continuous_vars  = NULL,
                categorical_vars = NULL,
                offset_vars      = NULL,
                C_hat            = NULL,
                DELTA            = 1e-12,
                sample_F         = (likelihood!= "Gaussian"),
                parallel         = FALSE,
                skip_postproc    = FALSE,
                t_test           = NULL,
                threshold        = 0.95,
                ...)
{
  # Create the model
  model <- lgp_model(formula       = formula,
                     data          = data,
                     likelihood    = likelihood,
                     prior         = prior,
                     uncertain_diagnosis = uncertain_diagnosis,
                     equal_effect  = equal_effect,
                     C_hat         = C_hat,
                     DELTA         = DELTA,
                     sample_F      = sample_F,
                     t_test        = t_test,
                     id_variable      = id_variable,
                     time_variable    = time_variable,
                     disAge_variable  = disAge_variable,
                     continuous_vars  = continuous_vars,
                     categorical_vars = categorical_vars,
                     offset_vars      = offset_vars)
  
  show(model)
  
  # Fit the model
  fit <- lgp_fit(model         = model, 
                 threshold     = threshold,
                 parallel      = parallel, 
                 skip_postproc = skip_postproc, 
                 ...)
  
  # Return the lgpfit object
  return(fit)
}
                      

#' Create an lgp model
#'
#' @export
#' @description Creates an object of class \code{lgpmodel}
#' @param formula A formula of the form
#' \code{y ~ } \code{x1 + x2 + x3}
#' defining the response variable \code{y} and covariates \code{xi}. All variables that
#' appear in the formula must exist as columns of \code{data}.
#' @param data A data frame containing (at least) the variables given in \code{formula}.
#' @param likelihood Either \code{"Gaussian"} (default), \code{"Poisson"} or \code{"NB"}.
#' @param prior Prior distribution. Can be created for example using the function 
#' \code{\link{prior_default}}.
#' @param uncertain_diagnosis Do we wish to model uncertainty in the disease onset?
#' @param equal_effect Is the disease effect assumed to be equally strong for all diseased individuals?
#' @param DELTA the amount of added jitter to ensure positive definiteness of the kernel
#' @param C_hat This can only be given if likelihood is not Gaussian. The signal \code{f} 
#' will the be transformed so that \code{g = } \code{exp(C_hat + f)}. If \code{NULL}, it will be
#' set to \code{C_hat = } \code{log(mean(y))}, where \code{y} is the response variable.
#' @param sample_F Determines if the function values are be sampled (must be \code{TRUE} if
#' likelihood is not Gaussian).
#' @param t_test Optional test time points. Should only be used if \code{sample_F = TRUE}.
#' Otherwise use \code{\link{lgp_predict}} after fitting the model.
#' @param id_variable Name of the unique subject identifier variable.
#' @param time_variable Name of the time variable.
#' @param disAge_variable Name of the disease-related age variable. If {NULL}, this
#' will be chosen to be "diseaseAge", if such covariate is found in the data.
#' @param continuous_vars Names of other continuous covariates. If \code{NULL}, the
#' remaining covariates that have floating point values are interpreted as continuous.
#' @param categorical_vars Names of categorical covariates that interact with the time variable.
#' If \code{NULL}, the remaining covariates that have integer values are interpreted 
#' as categorical.
#' @param offset_vars Names of the categorical covariates that are treated as 
#' time-independent group offsets. If \code{NULL}, no variables are interpreted as such
#' covariates.
#' @return An object of class \code{lgpmodel}.
#' @seealso For fitting the model, see \code{\link{lgp_fit}}.
lgp_model <- function(formula,
                      data,
                      likelihood       = "Gaussian",
                      prior            = prior_default(likelihood),
                      uncertain_diagnosis = FALSE,
                      equal_effect     = TRUE,
                      C_hat            = NULL,
                      DELTA            = 1e-12,
                      sample_F         = (likelihood!= "Gaussian"),
                      t_test           = NULL,
                      id_variable      = "id",
                      time_variable    = "age",
                      disAge_variable  = NULL,
                      continuous_vars  = NULL,
                      categorical_vars = NULL,
                      offset_vars      = NULL)
{
  # Model as a string
  fc <- as.character(formula)
  f  <- paste(fc[2], fc[1], fc[3])
  
  # Variable type info
  varInfo <- list(id_variable       = id_variable,
                  time_variable     = time_variable, 
                  disAge_variable   = disAge_variable, 
                  continuous_vars   = continuous_vars,
                  categorical_vars  = categorical_vars,
                  offset_vars       = offset_vars,
                  response_variable = fc[2])
  
  # Parse, check and preprocess the input
  PREPROC <- create_stan_input(formula        = formula,
                               data           = data,
                               prior          = prior,
                               likelihood     = likelihood,
                               varInfo        = varInfo,
                               standardize    = (likelihood %in% c("Gaussian", "none")),
                               uncertain_diagnosis = uncertain_diagnosis,
                               equal_effect   = equal_effect,
                               C_hat          = C_hat,
                               DELTA          = DELTA,
                               sample_F       = sample_F,
                               t_test         = t_test,
                               verbose        = TRUE)
  
  
  
  # Data to Stan
  stan_dat <- PREPROC$stan_dat
  
  # Create model info
  lh_str <- likelihood_as_str(stan_dat$LH)
  info <- list(likelihood      = lh_str,
               formula         = f,
               varInfo         = PREPROC$varInfo,
               sample_F        = as.logical(stan_dat$F_is_sampled),
               C_hat           = stan_dat$C_hat,
               DELTA           = stan_dat$DELTA,
               component_names = lgp_component_names(stan_dat),
               covariate_names = lgp_covariate_names(stan_dat),
               response_name   = fc[2])
  
  # Create the 'lgpmodel' object
  out <- new("lgpmodel",
             data     = data,
             stan_dat = stan_dat,
             scalings = PREPROC$scalings,
             info     = info
  )
  return(out)
}


#' Fit an lgp model
#'
#' @export
#' @description Samples the posterior of an additive Gaussian process regression model
#' using \code{\link{rstan}}.
#' @param model An object of class \code{lgpmodel}.
#' @param ... Optional arguments passed to \code{rstan::sampling}, for example \code{iter},
#' \code{chains} or \code{control}. See \code{\link[rstan]{sampling}} for the possible arguments.
#' @param parallel Determines if the chain will be run in parallel (\code{default = FALSE}).
#' If \code{TRUE}, then Stan is run by first defining 
#' \code{options(mc.cores = } \code{parallel::detectCores())}.
#' @param skip_postproc In this mode the postprocessing after running Stan is skipped.
#' @param threshold Covariate selection threshold.
#' @return An object of class \code{lgpfit}.
#' @seealso For the possible additional arguments, see \code{\link[rstan]{sampling}}.
#' For creating the \code{lgpmodel} input, see \code{\link{lgp_model}}.
#'
lgp_fit <- function(model, threshold, parallel = FALSE, skip_postproc = FALSE, ...)
{
  
  # Set Stan auto_write to true
  rstan::rstan_options(auto_write = TRUE)
  
  # Check if parallelization should be used
  if(parallel){
    options(mc.cores = parallel::detectCores())
  }else{
    options(mc.cores = 1)
  }
  
  # Run stan
  stan_dat  <- model@stan_dat
  stan_fit  <- rstan::sampling(object = stanmodels[["lgp"]], data = stan_dat, ...)
  
  # Initialize the 'lgpfit' object
  fit <- new("lgpfit", stan_fit = stan_fit, model = model)
  

  # Finalize the 'lgpfit' object
  tryCatch({
    if(!skip_postproc){
      fit@Rhat <- assess_convergence(fit)
      fit      <- postproc(fit, threshold, average_before_variance = F)
      show(fit)
    }else{
      cat("* Skipped postprocessing.\n")
    }
  }, error = function(e) {
    warning(e)
    cat("* Postprocessing failed, the error is printed as a warning.\n")
  }
  )
  
  # Return
  return(fit)
}


#' Compute predictions for a fitted model
#'
#' @export
#' @description Compute predictions for a fitted model. Only possible for models with
#' Gaussian likelihood.
#' @param fit An object of class \code{lgpfit}.
#' @param X_test The test points where the predictions should be computed.
#' @param samples The predictions can be computed either by using only the posterior mean
#' \cr (\code{samples="mean"}) or median (\code{samples="median"}) parameters, or by averaging over
#' all parameter samples (\code{samples="all"}). This can also be a set of indices, for example 
#' \code{samples=c(1:10)} gives the averaged predictions over the parameter samples 1...10.
#' @param print_progress Should progress be printed (if there is more than one sample)?
#' @param print_parameters Should the parameters be printed?
#' @return A data frame.
#' @seealso 
#' \itemize{
#' \item For creating an \code{lgpfit} object, see \code{\link{lgp_fit}}.
#' \item For creating an \code{lgpmodel} object, see \code{\link{lgp_model}}.
#' }
lgp_predict <- function(fit, 
                        X_test, 
                        samples = "mean",
                        print_progress = TRUE,
                        print_parameters = FALSE)
{
  
  # Check input correctness
  if(class(fit)!="lgpfit") stop("Class of 'fit' must be 'lgpfit'!")
  model <- fit@model
  if(class(model)!="lgpmodel") stop("Class of 'fit@model' must be 'lgpmodel'!")
  info     <- model@info
  cnames   <- info$component_names
  stan_dat <- model@stan_dat
  n        <- stan_dat$n
  X_sd     <- t(stan_dat$X)
  X_sd     <- X_sd[1:n,]
  Xnn_sd   <- t(stan_dat$X_notnan)
  Xnn_sd   <- Xnn_sd[1:n] 
  D        <- stan_dat$D
  LH       <- stan_dat$LH
  
  
  # Other model info
  fields   <- c("HMGNS", "UNCRT", "caseID_to_rows", "row_to_caseID", "DELTA")
  info     <- stan_dat[fields]
  
  if(info$UNCRT==1){
    stop("lgp_predict not implemented to work with UNCRT=1 yet!")
  }
  
  # Input checking
  if(LH!=1){
    stop("Computing predictions outside the data is possible only for models",
         " with Gaussian likelihood.")
  }
  cn1  <- colnames(X_test)
  cn2  <- colnames(X_sd)
  if(length(cn1) != length(cn2)){
    cn_str <- paste(colnames(X_sd), collapse=", ")
    stop("X_test must be a data frame with column names {", cn_str, "} (in this order).", sep="")
  }
  if(!all(cn1==cn2)){
    cn_str <- paste(colnames(X_sd), collapse=", ")
    stop("X_test must be a data frame with column names {", cn_str, "} (in this order).", sep="")
  }
  if(class(model)!="lgpmodel"){
    stop("model must be an object of class lgpmodel!")
  }
  if(class(fit)!="lgpfit"){
    stop("fit must be an object of class lgpfit!")
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
  y_data <- as.numeric(stan_dat$y)
  
  # All option to sample indicides
  if(is.character(samples)){
    if(samples=="all"){
      PAR <- hyperparam_samples(fit)
      samples <- 1:dim(PAR)[1]
    }
  }
  
  # Get parameters and compute predictions
  cat("* Computing predictions ")
  
  if(is.character(samples)){
    
    # Predictions using a single parameter estimate
    if(samples %in% c("mean", "median")){
      
      # Get parameter estimate
      params <- hyperparam_estimate(fit, samples)
      cat("using posterior ", samples, " parameters. \n", sep = "")
      if(print_parameters){print(params)}
      
      # Computation
      F_test <- compute_predictions(X_data, y_data, X_test, params, D, info, cnames)
      F_mu   <- F_test$F_mu
      F_var  <- F_test$F_var
      
    }else{
      stop("\ninvalid 'samples' input for lgp_predict!")
    }
    
  }else{
    # Averaged predictions using multiple samples
    cat("using ", length(samples), " parameter samples. \n\n", sep = "")
    PAR <- hyperparam_samples(fit, samples)
    
    # Print parameters
    if(print_parameters){print(PAR)}
    
    # Dimensions and sanity check
    ns  <- dim(PAR)[1]
    if(ns!=length(samples)){stop('sanity check failed!')}
    p <- dim(X_test)[1]
    d <- dim(X_test)[2]
    
    for(i_smp in 1:ns){
      
      # Get current parameter sample
      tmp <- colnames(PAR)
      params <- as.numeric(PAR[i_smp,])
      names(params) <- tmp
      
      # Computation
      F_test <- compute_predictions(X_data, y_data, X_test, params, D, info, cnames)
      
      # Averaging
      if(i_smp==1){
        F_mu   <- 1/ns * F_test$F_mu
        F_var  <- 1/ns * F_test$F_var
      }else{
        F_mu   <- F_mu + 1/ns * F_test$F_mu
        F_var  <- F_var + 1/ns * F_test$F_var
      }
      
      # Print  progress
      if(print_progress){
        cat('.')
        if(i_smp %% 10 == 0){ cat(' ')  }
        if(i_smp %% 50 == 0){ cat(i_smp, '\n') }
        if(i_smp %% 200 == 0){ cat('\n')}
        if(i_smp == ns){ cat('\n')}
      }
    }
    
  }
  
  # Return
  cat("* Done!\n\n")
  ret <- list(F_mu   = F_mu, 
              F_var  = F_var, 
              X_test_scaled = X_test,
              params = params)
  return(ret)
  
}

