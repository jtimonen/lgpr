#' The 'lgpr' package.
#'
#' @description
#' \itemize{
#'   \item A package for Bayesian additive Gaussian process (GP) modeling with
#' interpretable kernels designed for longitudinal data.
#'   \item Features inference of covariate effects
#'   and covariate relevance assesment.
#'   \item Models can include non-stationary, heterogeneous and temporally
#' uncertain effects.
#'   \item Bayesian inference of the model (hyper)parameters using
#' \code{\link[rstan]{rstan}}.
#'   \item  Functions for visualizing longitudinal data,
#' posterior draws, model predictions and inferred covariate effects are
#' also provided.
#' }
#'
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords Gaussian processes, longitudinal data, Stan, covariate relevances,
#' interpretable models
#'
#' @section Overview:
#' See the documentation of
#' \itemize{
#'  \item \code{\link{lgp}} for info on how to specify and fit models
#'  \item \code{\link{ppc}} for prior and posterior predictive checks
#'  \item \code{\link{relevances}} for component/covariate relevance assessment
#'  \item \code{\link{select}} for component/covariate selection
#'  \item \code{\link{get_pred}} for model predictions and inferred components
#'  at data points
#'  \item \code{\link{pred}} for computing out-of-sample predictions and
#'  inferred components
#'  \item \code{\link{plot_pred}} and \code{\link{plot_f}} for visualizing
#'  predictions and inferred components
#' }
#'
#' @section Tutorials:
#' See tutorials at \url{https://jtimonen.github.io/lgpr-usage/index.html}.
#'
#' @section Citation:
#' Run \code{citation("lgpr")} to get citation information.
#'
#' @section Feedback:
#' Bug reports, PRs, enchancement ideas or user experiences in general are
#' welcome and appreciated. Create an issue in Github or email the author.
#'
#' @docType package
#' @name lgpr-package
#' @aliases lgpr
#' @useDynLib lgpr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan get_rng get_stream
#'
#' @references
#' \enumerate{
#'   \item Carpenter, B. et al. (2017).
#'   \emph{Stan: A probabilistic programming language}. Journal of Statistical
#'    Software 76(1).
#'   \item Jonah Gabry, Ben Goodrich and Martin Lysy (2019).
#'   \emph{rstantools: Tools for Developing R Packages Interfacing with 'Stan'}.
#'   R package version 2.0.0.
#'   \item Gabry, J. and Mahr, T. (2019).
#'   \emph{bayesplot: Plotting for Bayesian Models}. R package version 1.7.0,
#'    http://mc-stan.org/bayesplot.
#'   \item Stan Development Team (2019). \emph{RStan: the R interface to Stan.}
#'   R package version 2.19.2. http://mc-stan.org/.
#' }
#'
NULL

#' Main function of the lgpr package
#'
#' @export
#' @description
#' Creates a model by calling \code{\link{create_model}} and fits by calling
#' \code{\link{sample_model}}.
#' @inheritParams create_model
#' @inheritParams sample_model
#' @inheritParams optimize_model
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
#' The \code{priors} argument must be a named list, like
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
                ...) {

  # Create model
  model <- create_model(
    formula, data, likelihood, prior, c_hat, num_trials, options,
    prior_only, verbose, sample_f
  )

  # Fit model
  if (verbose) cat("\nSampling model...\n")
  fit <- sample_model(model = model, ...)
  if (verbose) cat("Done.\n")
  return(fit)
}

#' Fitting a model
#'
#' @description
#' \itemize{
#'   \item The \code{sample_model} function takes an \linkS4class{lgpmodel}
#'   object, fits it using  \code{\link[rstan]{sampling}}, and
#'   returns an object of class \linkS4class{lgpfit}.
#'   \item The \code{optimize_model} function takes an \linkS4class{lgpmodel}
#'   object and fits it using \code{\link[rstan]{optimizing}}.
#' }
#' @name sample_model
#' @family main functions
NULL

#' @rdname sample_model
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param ... Optional arguments passed to
#' \code{\link[rstan]{sampling}} or \code{\link[rstan]{optimizing}}.
sample_model <- function(model, ...) {
  num_obs <- get_num_obs(model)
  large_data_msg(num_obs, 300)
  object <- stanmodels[[model@stan_model_name]]
  data <- model@stan_input
  sfit <- rstan::sampling(object = object, data = data, check_data = TRUE, ...)
  if (sfit@mode == 2) {
    print(sfit)
    stop("Failed to create fit.")
  }
  new("lgpfit", model = model, stan_fit = sfit)
}

#' @rdname sample_model
#' @export
optimize_model <- function(model, ...) {
  object <- stanmodels[[model@stan_model_name]]
  data <- model@stan_input
  rstan::optimizing(object = object, data = data, check_data = TRUE, ...)
}

#' Set the GP mean vector, taking TMM or other normalization
#' into account
#'
#' @export
#' @description Creates the \code{c_hat} input for \code{lgp},
#' so that it accounts for normalization between data points in the
#' \code{"poisson"} or \code{"nb"} observation model
#' @param y response variable, vector of length \code{n}
#' @param norm_factors normalization factors, vector of length \code{n}
#' @return a vector of length \code{n}, which can be used as
#' the \code{c_hat} input to the \code{lgp} function
adjusted_c_hat <- function(y, norm_factors) {
  check_lengths(y, norm_factors)
  check_non_negative_all(y)
  check_integer_all(y)
  check_positive_all(norm_factors)
  c_hat <- log(mean(y))
  c_hat <- rep(c_hat, length(y))
  c_hat <- c_hat + log(norm_factors)
  return(c_hat)
}
