#' Create input for Stan
#'
#' @description Parses the \code{formula} and \code{data} input to
#' \code{\link{lgp_model}}. Also performs many input checks.
#' @inheritParams lgp_model
#' @param standardize Should the response variable be standardized?
#' @param varInfo Variable type info.
#' @param verbose Should more verbose output be printed?
#' @return A list containing the data to be given to \code{rstan::sampling},
#' some info about preprocessing and all the information about scaling the
#' inputs and response, and updated variable type info.
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
                              verbose,
                              variance_mask,
                              N_trials,
                              skip_gen_quant) {

  # Parse likelihood
  LH <- likelihood_as_int(likelihood)

  # Check sample_F input
  if (!sample_F && !(LH %in% c(0, 1))) {
    stop("sample_F must be true if likelihood is not Gaussian or none")
  }

  # Check the model formula
  check_formula(formula, data)

  # Get the (scaled) response variable
  RESP <- get_response(data, varInfo, standardize, LH)
  response <- RESP$response

  # Resolve covariate types
  TYP <- check_data(data, varInfo, verbose)
  types <- TYP$types
  varInfo <- TYP$varInfo

  # Take only those covariates that are needed by Stan
  covariates <- create_covariates_stan(data, varInfo, types, formula, verbose)
  D <- covariates$D
  X <- covariates$X
  X_notnan <- covariates$X_notnan

  # Compute dimensions and standardize covariates
  stan_dims <- get_model_dims(X, D)
  N_cat <- set_N_cat(X, D)
  SCL <- standardize_inputs(X, D)
  X <- SCL$X

  # Check and possibly edit N_trials, C_hat and norm_factors
  N_trials <- set_N_trials(N_trials, response, LH)
  C_hat <- set_C_hat(C_hat, response, LH, N_trials)

  # Check that variable types make sense
  check_varInfo(varInfo)

  # Create the list that is the Stan input
  stan_dat <- list(
    X = t(X),
    X_notnan = X_notnan,
    y = response,
    y_int = round(response),
    LH = LH,
    N_trials = N_trials,
    N_cat = N_cat,
    C_hat = C_hat,
    F_IS_SAMPLED = as.numeric(sample_F),
    USE_VAR_MASK = as.numeric(variance_mask),
    VERBOSE = as.numeric(verbose),
    DELTA = DELTA,
    SKIP_GQ = as.numeric(skip_gen_quant)
  )

  # Get some variables related to diseased individuals
  stan_dis <- get_diseased_info(
    D, X, X_notnan,
    uncertain_effect_time,
    equal_effect,
    SCL$TSCL
  )

  # Parse the prior
  stan_prior <- prior_to_stan(
    D, prior,
    stan_dis$HMGNS,
    stan_dis$UNCRT,
    stan_dis$N_cases,
    stan_dis$T_observed,
    stan_dis$T_last
  )

  # Return
  ret <- list(
    stan_dat = c(stan_dat, stan_dims, stan_dis, stan_prior),
    varInfo = varInfo,
    scalings = list(TSCL = SCL$TSCL, CSCL = SCL$CSCL, YSCL = RESP$SCL)
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
create_covariates_stan <- function(data, varInfo, types, formula, verbose) {

  # Create the design matrix X
  XD <- stan_input_X_and_D(data, varInfo, types, formula, verbose)

  # Do this because as.matrix() sometimes creates character matrices
  X <- data.matrix(XD$X)
  D <- XD$D

  # Get location of NaNs because Stan does not accept NaNs
  if (D[3] > 0) {
    X_notnan <- 1 - is.nan(X[, 3])
    X[which(X_notnan == 0), 3] <- 0
  } else {
    X_notnan <- rep(0, length(X[, 1]))
  }
  X_notnan <- as.vector(X_notnan)

  # If there are still NaNs, throw an error
  n_nans <- sum(is.nan(X))
  if (n_nans > 0) {
    stop("Only the diseaseAge column of the data can contain NaNs!")
  }

  # Return
  covariates <- list(
    X = X,
    X_notnan = X_notnan,
    D = D
  )

  return(covariates)
}
