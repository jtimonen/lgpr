
#' Create an lgp model
#'
#' @export
#' @description Creates an object of class \code{lgpmodel}
#' @param formula A formula specified using the common \code{\link{formula}}
#' syntax, such as \code{y ~ x1 + x2:z1 + } \code{x2:z2 + z2}.
#' \itemize{
#'   \item The formula must contain exatly one tilde (\code{~}), with response
#'   variable on the left-hand side and model terms on the right-hand side.
#'   \item Terms are be separated by a plus (\code{+}) sign.
#'   \item Terms can consits of a single variable name, such as \code{x}, or
#'   an interaction of two variables, such as \code{x:z}.
#'   \item In single-variable terms, the variable can be either continuous or
#'   categorical, whereas in interaction terms the variable
#'   on the left-hand side of the colon (\code{:}) has to be continuous and the
#'   one on the right-hand side has to be categorical (a factor).
#'   \item All variables appearing in \code{formula} must be
#'   found in \code{data}.
#' }
#' @param data A data frame where each column corresponds to one variable, and
#' each row is one observation. Continuous covariates must have type
#' \code{"numeric"} and categorical ones must have type \code{"factor"}.
#' Note that you need to indicate missing values of \code{disAge_variable} with
#' \code{NaN}.
#' @param likelihood Determines the observation model. Must be either
#' \code{"Gaussian"} (default), \code{"Poisson"}, \code{"NB"} (negative
#' binomial) or \code{"binomial"}. To
#' use Bernoulli likelihood, use \code{likelihood="binomial"} and set
#' \code{N_trials} as a vector of ones.
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters. It is recommended to first create this using the function
#' \code{\link{prior_default}}, and then possibly modify it.
#' @param uncertain_effect_time Do we wish to model uncertainty in the disease
#' effect time?
#' @param equal_effect Is the disease effect assumed to be equally strong for
#' all diseased individuals?
#' @param DELTA the amount of added jitter to ensure positive definiteness of
#' the kernel
#' @param C_hat The GP mean. Must be a vector of length \code{dim(data)[1]}, or
#' a real number defining a constant GP mean. If \code{NULL}, this is set to
#'  \itemize{
#'    \item \code{C_hat = 0}, if \code{likelihood} is \code{"Gaussian"}, because
#'    with
#'    Gaussian likelihood the response variable is by default centered to have
#'    zero mean.
#'    \item \code{C_hat = } \code{log(mean(y))} if \code{likelihood} is
#'    \code{"Poisson"} or \code{"NB"},
#'    \item \code{C_hat = } \code{log(p/(1-p))}, where
#'    \code{p = mean(y/N_trials)} if \code{likelihood} is \code{"binomial"},
#'  }
#' where \code{y} denotes the response variable. You can modify this vector to
#' account for normalization between data points. With Gaussian likelihood
#' though, do not modify this argument, normalize the data beforehand instead.
#' @param sample_F Determines if the function values are be sampled (must be
#' \code{TRUE} if likelihood is not \code{"Gaussian"}).
#' @param id_variable Name of the unique subject identifier variable (default =
#' \code{"id"}).
#' @param time_variable Name of the time variable (default = \code{"age"}).
#' @param disAge_variable Name of the disease-related age variable. If
#' \code{NULL}, this will be chosen to be \code{"diseaseAge"}, if such covariate
#' is found in the data. You need to indicate missing values of
#' \code{disAge_variable} with \code{NaN} in \code{data}.
#' @param offset_vars Names of the categorical covariates that are treated as
#' time-independent group offsets. If \code{NULL} (default), no variables are
#' interpreted as such covariates.
#' @param variance_mask Should a variance mask be used to force disease
#' component variance to zero before disease onset?
#' @param N_trials This argument (number of trials) is only needed when
#' likelihood is binomial. Must have length one or equal to number of data
#' points. Setting \code{N_trials=1} corresponds to Bernoulli observation model.
#' @param skip_gen_quant If this is true, the generated quantities block of Stan
#' is skipped.
#' @param verbose Should more verbose output be printed?
#' @return An object of class \code{lgpmodel}.
#' @seealso For fitting the model, see \code{\link{lgp_fit}}.
lgp_model <- function(formula,
                      data,
                      likelihood = "Gaussian",
                      prior = prior_default(),
                      uncertain_effect_time = FALSE,
                      equal_effect = TRUE,
                      C_hat = NULL,
                      DELTA = 1e-8,
                      sample_F = NULL,
                      id_variable = "id",
                      time_variable = "age",
                      disAge_variable = NULL,
                      continuous_vars = NULL,
                      categorical_vars = NULL,
                      offset_vars = NULL,
                      variance_mask = TRUE,
                      N_trials = NULL,
                      skip_gen_quant = FALSE,
                      verbose = FALSE) {
  # Model as a string
  fc <- as.character(formula)
  f <- paste(fc[2], fc[1], fc[3])

  # Is F sampled
  likelihood <- tolower(likelihood)
  lh_gauss_or_none <- likelihood %in% c("gaussian", "none")
  if (is.null(sample_F)) {
    sample_F <- !lh_gauss_or_none
  }

  # Variable type info
  varInfo <- list(
    id_variable = id_variable,
    time_variable = time_variable,
    disAge_variable = disAge_variable,
    continuous_vars = continuous_vars,
    categorical_vars = categorical_vars,
    offset_vars = offset_vars,
    response_variable = fc[2]
  )

  # Parse, check and preprocess the input
  PREPROC <- create_stan_input(
    formula = formula,
    data = data,
    prior = prior,
    likelihood = likelihood,
    varInfo = varInfo,
    standardize = lh_gauss_or_none,
    uncertain_effect_time = uncertain_effect_time,
    equal_effect = equal_effect,
    C_hat = C_hat,
    DELTA = DELTA,
    sample_F = sample_F,
    verbose = verbose,
    variance_mask = variance_mask,
    N_trials = N_trials,
    skip_gen_quant = skip_gen_quant
  )

  # Data to Stan
  stan_dat <- PREPROC$stan_dat

  # Create model info
  lh_str <- likelihood_as_str(stan_dat$LH)
  info <- list(
    likelihood = lh_str,
    formula = f,
    varInfo = PREPROC$varInfo,
    sample_F = as.logical(stan_dat$F_IS_SAMPLED),
    C_hat = stan_dat$C_hat,
    DELTA = stan_dat$DELTA,
    component_names = lgp_component_names(stan_dat),
    covariate_names = lgp_covariate_names(stan_dat),
    response_name = fc[2],
    variance_mask = variance_mask,
    N_trials = N_trials
  )

  # Create the 'lgpmodel' object
  out <- new("lgpmodel",
    data = data,
    stan_dat = stan_dat,
    scalings = PREPROC$scalings,
    info = info
  )
  return(out)
}
