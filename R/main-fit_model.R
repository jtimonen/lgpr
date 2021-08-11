#' Main function of the 'lgpr' package
#'
#' @export
#' @description
#' Creates an additive Gaussian process model using
#' \code{\link{create_model}} and fits it using \code{\link{sample_model}}.
#' See the
#' \href{https://jtimonen.github.io/lgpr-usage/articles/math.html}{Mathematical description of lgpr models}
#' vignette for more information about the connection between different options
#' and the created statistical model.
#'
#' @inheritParams create_model
#' @inheritParams sample_model
#' @inheritParams optimize_model
#' @param verbose Can messages be printed during model creation? Has no
#' effect if \code{quiet=TRUE}.
#'
#' @section Model formula syntax:
#' There are two ways to define the model formula:
#' \enumerate{
#'   \item Using a common \code{\link[stats]{formula}}-like syntax, like in
#'   \code{y ~ age +} \code{age|id} \code{ + sex}. Terms can consist of a
#'   single variable, such as \code{age}, or an interaction of two variables,
#'   such as \code{age|id}. In single-variable terms, the variable can be either
#'   continuous (numeric) or categorical (factor), whereas in interaction terms
#'   the variable on the left-hand side of the vertical bar (\code{|}) has to
#'   be continuous and the one on the right-hand side has to be categorical.
#'   Formulae specified using this syntax are translated to the advanced format
#'   so that
#'   \itemize{
#'     \item single-variable terms become \code{gp(x)} if
#'     variable \code{x} is numeric and \code{zs(x)} if \code{x} is a factor
#'     \item interaction terms \code{x|z} become \code{gp(x)*zs(z)}
#'   }
#'   \item Using the advanced syntax, like in \code{y ~ gp(age) +}
#'   \code{gp(age)*zs(id) +} \code{het(id)*gp_vm(disAge)}.
#'   This creates \linkS4class{lgprhs} objects, which consist of
#'  \linkS4class{lgpterm}s, which consist of \linkS4class{lgpexpr}s.
#'  This approach must be used if creating nonstationary, heterogeneous or
#'  temporally uncertain components.
#' }
#' Either one of the approaches should be used and they should not be mixed.
#'
#' @section Defining priors:
#' The \code{prior} argument must be a named list, like
#' \code{list(alpha=student_t(4), wrp=igam(30,10))}. See examples in tutorials.
#' Possible allowed names are
#' \itemize{
#'  \item \code{"alpha"} = component magnitude parameters
#'  \item \code{"ell"} = component lengthscale parameters
#'  \item \code{"wrp"} = input warping steepness parameters
#'  \item \code{"sigma"} = noise magnitude (Gaussian obs. model)
#'  \item \code{"phi"} = inv. overdispersion (negative binomial obs. model)
#'  \item \code{"gamma"} = overdispersion (beta-binomial obs. model)
#'  \item \code{"beta"} = heterogeneity parameters
#'  \item \code{"effect_time"} = uncertain effect time parameters
#'  \item \code{"effect_time_info"} = additional options for the above
#' }
#' See \code{\link{priors}} for functions that can be
#' used to define the list elements. If a parameter of a model is not given
#' in this list, a default prior will be used for it.
#'
#' @section When to not use default priors:
#'
#' It is not recommended to use default priors blindly. Rather, priors should
#' be specified according to the knowledge about the problem at hand, as in any
#' Bayesian analysis. In \code{lgpr} this is especially important when
#' \enumerate{
#'  \item Using a non-Gaussian likelihood or otherwise setting
#'  \code{sample_f = TRUE}. In this case the response variable is not
#'  normalized, so the scale on which the data varies must be taken into
#'  account when defining priors of the signal magnitude parameters
#'  \code{alpha} and possible noise parameters (\code{sigma}, \code{phi},
#'  \code{gamma}). Also it should be checked if \code{c_hat} is set in a
#'  sensible way.
#'  \item Using a model that contains a \code{gp_ns(x)} or \code{gp_vm(x)}
#'  expression in its formula. In this case the corresponding covariate
#'  \code{x} is not normalized, and the prior for the input warping steepness
#'  parameter \code{wrp} must be set according to the expected width of the
#'  window in which the nonstationary effect of \code{x} occurs. By default,
#'  the width of this window is about 36, which has been set assuming that
#'  the unit of \code{x} is months.
#' }
#'
#'
#' @name lgp
#' @family main functions
#' @return Returns an object of the S4 class \linkS4class{lgpfit}.
NULL

#' @rdname lgp
lgp <- function(formula,
                data,
                likelihood = "gaussian",
                prior = NULL,
                c_hat = NULL,
                num_trials = NULL,
                options = NULL,
                prior_only = FALSE,
                verbose = FALSE,
                sample_f = !(likelihood == "gaussian"),
                quiet = FALSE,
                skip_postproc = sample_f,
                ...) {
  if (quiet) verbose <- FALSE

  # Create model
  log_progress("Creating model...", verbose)
  model <- create_model(
    formula, data, likelihood, prior, c_hat, num_trials, options,
    prior_only, verbose, sample_f
  )
  if (verbose) {
    log_progress("\nModel created, printing it here.")
    print(model)
  }

  # Fit model
  log_progress("\nSampling model...", verbose)
  fit <- sample_model(
    model = model,
    verbose = verbose,
    quiet = quiet,
    skip_postproc = skip_postproc,
    ...
  )
  return(fit)
}

#' Fitting a model
#'
#' @description
#' \itemize{
#'   \item \code{sample_model} takes an \linkS4class{lgpmodel}
#'   object and fits it using \code{\link[rstan]{sampling}}.
#'   \item \code{optimize_model} takes an \linkS4class{lgpmodel}
#'   object and fits it using \code{\link[rstan]{optimizing}}.
#' }
#' @name sample_model
#' @family main functions
#' @return
#' \itemize{
#'   \item \code{sample_model} returns an object of class \linkS4class{lgpfit}
#'   containing the parameter draws, the original \code{model} object,
#'   and possible postprocessing results. See documentation of
#'   \linkS4class{lgpfit} for more information.
#'   \item \code{optimize_model} directly returns the list returned by
#' \code{\link[rstan]{optimizing}}. See its documentation for more information.
#' }
NULL

#' @rdname sample_model
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param quiet Should all output messages be suppressed? You need to set
#' also \code{refresh=0} if you want to suppress also the progress update
#' messages from \code{\link[rstan]{sampling}}.
#' @param skip_postproc Should all postprocessing be skipped? If this is
#' \code{TRUE}, the returned \linkS4class{lgpfit} object will likely be
#' much smaller (if \code{sample_f=FALSE}).
#' @param verbose Can messages be printed?
#' @param ... Optional arguments passed to
#' \code{\link[rstan]{sampling}} or \code{\link[rstan]{optimizing}}.
sample_model <- function(model, verbose = TRUE, quiet = FALSE,
                         skip_postproc = is_f_sampled(model), ...) {
  if (quiet) verbose <- FALSE
  num_obs <- get_num_obs(model)
  LIMIT <- 300 # __HARDCODED__
  large_data_msg(num_obs, LIMIT)
  object <- get_stan_model(model)
  data <- model@stan_input

  # Run sampling
  stan_fit <- rstan::sampling(
    object = object,
    data = data,
    check_data = TRUE,
    pars = "eta",
    include = FALSE,
    ...
  )
  log_progress("Sampling done.", verbose)
  if (stan_fit@mode == 2) {
    if (!quiet) print(stan_fit)
    stop("Failed to create stanfit.")
  }

  # Create the lgpfit object
  fit <- new("lgpfit",
    model = model,
    stan_fit = stan_fit,
    num_draws = get_num_postwarmup_draws(stan_fit),
    postproc_results = list()
  )

  # Postprocess
  log_progress("\nPostprocessing...", verbose)
  if (skip_postproc) {
    return(fit)
  } else {
    verbose_postproc <- if (quiet) FALSE else TRUE
    tryCatch(
      {
        fit <- postproc(fit, verbose = verbose_postproc)
      },
      error = function(e) {
        warning("\nPostprocessing failed. Reason:")
        warning(e)
      }
    )
  }
  log_progress("Postprocessing done.", verbose)
  return(fit)
}


#' @rdname sample_model
#' @export
optimize_model <- function(model, ...) {
  num_obs <- get_num_obs(model)
  large_data_msg(num_obs, 300)
  object <- get_stan_model(model)
  data <- model@stan_input
  rstan::optimizing(object = object, data = data, check_data = TRUE, ...)
}
