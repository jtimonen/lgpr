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
                uncertain_effect_time = FALSE,
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
                variance_mask    = TRUE,
                ell_smooth       = "ell_shared",
                ell_smooth_multip = 1,
                cat_interact_kernel_type = "categorical",
                N_trials         = NULL,
                ...)
{
  # Create the model
  model <- lgp_model(formula          = formula,
                     data             = data,
                     likelihood       = likelihood,
                     prior            = prior,
                     uncertain_effect_time = uncertain_effect_time,
                     equal_effect     = equal_effect,
                     C_hat            = C_hat,
                     DELTA            = DELTA,
                     sample_F         = sample_F,
                     t_test           = t_test,
                     id_variable      = id_variable,
                     time_variable    = time_variable,
                     disAge_variable  = disAge_variable,
                     continuous_vars  = continuous_vars,
                     categorical_vars = categorical_vars,
                     offset_vars      = offset_vars,
                     variance_mask    = variance_mask,
                     cat_interact_kernel_type  = cat_interact_kernel_type,
                     N_trials         = N_trials)
  
  show(model)
  
  # Fit the model
  fit <- lgp_fit(model         = model, 
                 threshold     = threshold,
                 parallel      = parallel, 
                 skip_postproc = skip_postproc,
                 ell_smooth    = ell_smooth,
                 ell_smooth_multip = ell_smooth_multip,
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
#' @param likelihood Determines the observation model. Must be either \code{"Gaussian"} 
#' (default), \code{"Poisson"}, \code{"NB" (negative binomial) or "binomial"}.
#' @param prior Prior distribution. Can be created for example using the function 
#' \code{\link{prior_default}}.
#' @param uncertain_effect_time Do we wish to model uncertainty in the disease effect time?
#' @param equal_effect Is the disease effect assumed to be equally strong for all diseased individuals?
#' @param DELTA the amount of added jitter to ensure positive definiteness of the kernel
#' @param C_hat This can only be given if likelihood is Poisson or NB. The signal \code{f} 
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
#' @param variance_mask Should a variance mask be used to force disease component
#' variance to zero before disease onset?
#' @param cat_interact_kernel_type Kernel type for categorical variables (other than id). 
#' Possible options are \code{"categorical"} (default) and \code{"binary"} (mask kernel where
#' only category \code{"1"} will have an effect).
#' @param N_trials This argument (number of trials) is only needed when likelihood is binomial.
#' Must have length one or equal to number of data points. Setting \code{N_trials=1} corresponds to 
#' Bernoulli observation model.
#' @return An object of class \code{lgpmodel}.
#' @seealso For fitting the model, see \code{\link{lgp_fit}}.
lgp_model <- function(formula,
                      data,
                      likelihood       = "Gaussian",
                      prior            = prior_default(likelihood),
                      uncertain_effect_time = FALSE,
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
                      offset_vars      = NULL,
                      variance_mask    = TRUE,
                      cat_interact_kernel_type = "categorical",
                      N_trials = NULL
)
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
                               uncertain_effect_time = uncertain_effect_time,
                               equal_effect   = equal_effect,
                               C_hat          = C_hat,
                               DELTA          = DELTA,
                               sample_F       = sample_F,
                               t_test         = t_test,
                               verbose        = TRUE,
                               variance_mask  = variance_mask,
                               cat_interact_kernel_type = cat_interact_kernel_type,
                               N_trials       = N_trials)
  
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
               response_name   = fc[2],
               variance_mask   = variance_mask,
               cat_interact_kernel_type = cat_interact_kernel_type,
               N_trials        = N_trials)
  
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
#' @param ell_smooth Defines how to determine smoothing lengthscale for corrected shared age 
#' effect inference. Possible options are
#' \enumerate{
#'    \item \code{"ell_shared"} (default) - the sampled lengthscale of the shared age 
#'    component is used as \code{ell_smooth}
#'    \item \code{"none"} - no correction will be performed
#'    \item A numeric argument that directly defines \code{ell_smooth}
#' }
#' @param ell_smooth_multip a multiplier for ell_smooth
#' @return An object of class \code{lgpfit}.
#' @seealso For the possible additional arguments, see \code{\link[rstan]{sampling}}.
#' For creating the \code{lgpmodel} input, see \code{\link{lgp_model}}.
#'
lgp_fit <- function(model, 
                    threshold, 
                    parallel          = FALSE, 
                    skip_postproc     = FALSE, 
                    ell_smooth        = "ell_shared",
                    ell_smooth_multip = 1,
                    ...)
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
  ver <- "NOT FOUND"
  tryCatch({
    ver <- get_pkg_description()$Version
  }, error = function(e){
    cat("unable to get version number\n")
  })
  fit <- new("lgpfit", stan_fit = stan_fit, model = model, 
             pkg_version = ver)
  

  # Finalize the 'lgpfit' object
  tryCatch({
    if(!skip_postproc){
      cat("* Begin postprocessing. \n")
      fit@Rhat <- assess_convergence(fit)
      fit      <- postproc(fit, 
                           threshold               = threshold, 
                           ell_smooth              = ell_smooth,
                           ell_smooth_multip       = ell_smooth_multip,
                           average_before_variance = FALSE)
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
#' \cr (\code{samples="mean"}), median (\code{samples="median"}), or MAP (\code{samples="map"}) 
#' parameters, or for all parameter samples (\code{samples="all"}). This can also be a set of indices, 
#' for example \code{samples=c(1:10)} gives predictions for the parameter samples 1...10.
#' @param print_progress Should progress be printed (if there is more than one sample)?
#' @param print_params Should the parameter values be printed? (only works if \code{samples}
#' is mean or median.)
#' @return A list.
#' @seealso 
#' \itemize{
#' \item For creating an \code{lgpfit} object, see \code{\link{lgp_fit}}.
#' \item For creating an \code{lgpmodel} object, see \code{\link{lgp_model}}.
#' }
lgp_predict <- function(fit, 
                        X_test, 
                        samples        = "mean",
                        print_progress = TRUE,
                        print_params   = FALSE)
{
  
  # Run some input checks and scale data correctly
  PP      <- predict_preproc(fit, X_test, samples)
  X_data  <- PP$X_data
  X_test  <- PP$X_test
  y_data  <- PP$y_data
  samples <- PP$samples
  info    <- PP$info
  D       <- PP$D
  cnames  <- PP$cnames
  LIST    <- list()
  TSCL    <- fit@model@scalings$TSCL
  
  cat("* Computing predictions ")
  
  if(is.character(samples)){
    
    # Predictions using a single parameter estimate
    if(samples %in% c("mean", "median", "map")){
      if(samples %in% c("mean", "median")){
        samples_str <- paste("posterior", samples)
      }else{
        samples_str <- "MAP"
      }
      cat("using ", samples_str, " parameters. \n", sep = "")
      params    <- hyperparam_estimate(fit, samples)
      if(print_params){
        print(params)
      }
      LIST[[1]] <- compute_predictions(X_data, y_data, X_test, params, D, info, cnames, TSCL)
    }else{
      stop("\ninvalid 'samples' input for lgp_predict! (", samples, ")")
    }
  }else{
    
    # Predictions using multiple samples
    cat("using ", length(samples), " parameter sample(s). \n\n", sep = "")
    PAR    <- hyperparam_samples(fit, samples)
    pnames <- colnames(PAR)
    ns     <- dim(PAR)[1]
    
    # Progress bar top
    if(print_progress){
      str    <- paste("|   ", seq(10,100,by=10), "%", sep="")
      top    <- paste(formatC(str, width = 4), collapse = " ")
      top    <- paste(top, "|")
      barlen <- nchar(top) - 1
      iprint <- ceiling(seq(1, ns, length.out = barlen))
      if(ns >=100){
        cat(top, "\n ")
      }
    }
    
    for(i_smp in 1:ns){
      # Predict with current parameter sample
      params        <- as.numeric(PAR[i_smp,])
      names(params) <- pnames
      LIST[[i_smp]] <- compute_predictions(X_data, y_data, X_test, params, D, info, cnames, TSCL)
      
      # Print progress
      if(print_progress %% ns > 1){
        if(i_smp %in% iprint){cat('=')}
        if(i_smp == ns){ cat('\n\n')}
      }
    }
    
  }
  
  # Return
  ret <- list(LIST = LIST, X_test_scaled = X_test)
  return(ret)
}


#' Compute predictions and log-posterior predictive density at test points
#'
#' @export
#' @description This is a convenience function that wraps \code{\link{lgp_predict}},
#' \code{\link{compute_lppd}} and \code{\link{plot_posterior_y}}.
#' @param fit an object of class \code{lgpfit}
#' @param test_data a test data matrix
#' @param verbose Should this print progress?
#' @param samples Sample indices or a keyword "mean", "median", "map", or "all".
#' @param plot should this return also a plot of the data and predictions?
#' @return a ggplot object or lppd
lgp_test <- function(fit, test_data, plot = FALSE, verbose = TRUE, samples = "mean"){
  
  # predict only at test points
  info    <- fit@model@info
  yname   <- info$response_name
  xnames  <- info$covariate_names
  idvar   <- info$varInfo$id_variable
  tvar    <- info$varInfo$time_variable
  xnames  <- union(c(idvar, tvar), xnames)
  dat     <- data.frame(test_data)
  X_test  <- dat[,xnames]
  y_test  <- dat[,yname]
  PRED    <- lgp_predict(fit, X_test, samples = samples, print_progress = verbose)
  LPPD    <- compute_lppd(PRED, y_test)
  mlppd   <- mean(as.numeric(LPPD))
  ret     <- list(lppd = LPPD, mlppd = mlppd)
  if(plot){
    p    <- plot_posterior_y(fit, PRED, test_data = test_data, uncertainty = "errorbar")
    subt <- paste("Mean log-posterior predictive density:", round(mlppd, 5))
    p    <- p + ggplot2::ggtitle(label = "Predictive distribution of y at test points",
                              subtitle = subt)
    ret$plot <- p
  }
  return(ret)
}
