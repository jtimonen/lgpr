# MAIN DOCUMENTATION PAGE -------------------------------------------------

#' The 'lgpr' package.
#'
#' @description Interpretable nonparametric modeling of longitudinal data
#' using additive Gaussian process regression. Contains functionality
#' for inferring covariate effects and assessing covariate relevances.
#' Models are specified using a convenient formula syntax, and can include
#' shared, group-specific, non-stationary, heterogeneous and temporally
#' uncertain effects. Bayesian inference for model parameters is performed
#' using 'Stan' (\code{\link[rstan]{rstan}}). The modeling approach and methods
#' are described in detail in
#' \href{https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab021/6104850}{Timonen et al. (2021)}.
#'
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords longitudinal additive model GP Stan interpretable Bayesian
#' covariate relevance
#'
#' @section Core functions:
#' Main functionality of the package consists of creating and fitting an
#' additive GP model:
#' \itemize{
#'  \item \code{\link{lgp}}: Specify and fit an additive GP model with one
#'  command.
#'  \item \code{\link{create_model}}: Define an \linkS4class{lgpmodel} object.
#'  \item \code{\link{sample_model}}: Fit a model by sampling the posterior
#'  distribution of its parameters and create an \linkS4class{lgpfit} object.
#'  \item \code{\link{pred}}: Computing model predictions and inferred
#'  covariate effects after fitting a model.
#'  \item \code{\link{relevances}}: Assessing covariate relevances after
#'  fitting a model.
#'  \item \code{\link{prior_pred}}: Prior predictive sampling to check
#'  if your prior makes sense.
#' }
#'
#' @section Visualization:
#' \itemize{
#'  \item \code{\link{plot_pred}}: Plot model predictions.
#'  \item \code{\link{plot_components}}: Visualize inferred model components.
#'  \item \code{\link{plot_draws}}: Visualize parameter draws.
#'  \item \code{\link{plot_data}}: Visualize longitudinal data.
#' }
#' @section Data:
#' The data that you wish to analyze with 'lgpr' should be in an \R
#' \code{data.frame} where columns correspond to measured variables and rows
#' correspond to observations. Some functions that can help working with such
#' data frames are:
#' \itemize{
#'  \item \code{\link{new_x}}: Creating new test points where the posterior
#'  distribution of any function component or sum of all components, or the
#'  posterior predictive distribution can be computed after model fitting.
#'  \item Other functions: \code{\link{add_factor}},
#'  \code{\link{add_factor_crossing}}, \code{\link{add_dis_age}},
#'  \code{\link{adjusted_c_hat}}.
#' }
#'
#' @section Vignettes and tutorials:
#' See \url{https://jtimonen.github.io/lgpr-usage/index.html}. The
#' tutorials focus on code and use cases, whereas the
#' \href{https://jtimonen.github.io/lgpr-usage/articles/math.html}{Mathematical description of lgpr models}
#' vignette describes the statistical models and how they can be customized in 'lgpr'.
#'
#' @section Citation:
#' Run \code{citation("lgpr")} to get citation information.
#'
#' @section Feedback:
#' Bug reports, PRs, enhancement ideas or user experiences in general are
#' welcome and appreciated. Create an issue in Github or email the author.
#'
#' @docType package
#' @name lgpr-package
#' @aliases lgpr
#' @useDynLib lgpr, .registration = TRUE
#' @import methods rstantools
#' @importFrom rstan get_stream
#' @importClassesFrom rstan stanfit
#'
#' @references
#' \enumerate{
#'   \item Timonen, J. et al. (2021).
#'   \emph{lgpr: an interpretable non-parametric method for inferring covariate
#'   effects from longitudinal data}. Bioinformatics,
#'   \href{https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab021/6104850}{url}.
#'   \item Carpenter, B. et al. (2017).
#'   \emph{Stan: A probabilistic programming language}. Journal of Statistical
#'    Software 76(1).
#' }
#'
NULL


# DATASETS ----------------------------------------------------------------

#' A very small artificial test data, used mostly for unit tests
#'
#' @format A data frame with 24 rows and 6 variables:
#' \describe{
#'   \item{id}{individual id, a factor with levels: 1, 2, 3, 4}
#'   \item{age}{age}
#'   \item{dis_age}{disease-related age}
#'   \item{blood}{a continuous variable}
#'   \item{sex}{a factor with 2 levels: Male, Female}
#'   \item{y}{a continuous variable}
#' }
#' @family built-in datasets
"testdata_001"

#' Medium-size artificial test data, used mostly for tutorials
#'
#' @format A data frame with 96 rows and 6 variables:
#' \describe{
#'   \item{id}{individual id, a factor with levels: 01-12}
#'   \item{age}{age}
#'   \item{diseaseAge}{disease-related age}
#'   \item{sex}{a factor with 2 levels: Male, Female}
#'   \item{group}{a factor with 2 levels: Case, Control}
#'   \item{y}{a continuous variable}
#' }
#' @family built-in datasets
#' @seealso \code{\link{read_proteomics_data}}
"testdata_002"


# S4 CLASSES --------------------------------------------------------------

#' Validate S4 class objects
#'
#' @param object an object to validate
#' @return \code{TRUE} if valid, otherwise reasons for invalidity
#' @name validate
NULL

#' @rdname validate
validate_lgpexpr <- function(object) {
  errors <- character()
  v0 <- nchar(object@covariate) > 0
  valid_funs <- c(
    "gp", "gp_ns", "gp_vm", "categ", "zs", "het", "unc"
  )
  v1 <- object@fun %in% valid_funs
  if (!v0) {
    errors <- c(errors, "covariate name cannot be empty")
  }
  if (!v1) {
    str <- paste0(valid_funs, collapse = ", ")
    msg <- paste0(
      "<fun> must be one of {", str,
      "}, found = '", object@fun, "'"
    )
    errors <- c(errors, msg)
  }
  return_true_or_errors(errors)
}

#' @rdname validate
validate_lgpformula <- function(object) {
  covs <- rhs_variables(object@terms)
  r <- object@y_name
  errors <- character()
  if (r %in% covs) {
    msg <- "the response variable cannot be also a covariate"
    errors <- c(errors, msg)
  }
  return_true_or_errors(errors)
}

#' @rdname validate
validate_lgpscaling <- function(object) {
  loc <- object@loc
  scale <- object@scale
  L1 <- length(loc)
  L2 <- length(scale)
  errors <- character()
  if (L1 != 1) {
    err <- paste("length(loc) should be 1, found = ", L1)
    errors <- c(errors, err)
  }
  if (L2 != 1) {
    err <- paste("length(scale) should be 1, found = ", L2)
    errors <- c(errors, err)
  }
  if (nchar(object@var_name) < 1) {
    err <- "variable name length must be at least 1 character"
    errors <- c(errors, err)
  }
  if (scale <= 0) {
    err <- paste0("scale must be positive, found = ", scale)
    errors <- c(errors, err)
  }
  return_true_or_errors(errors)
}

#' @rdname validate
validate_lgpfit <- function(object) {
  errors <- character()

  # Validate num_draws
  N1 <- nrow(get_draws(object, NULL, NULL, pars = c("alpha")))
  N2 <- object@num_draws
  msg <- paste0(
    "invalid num_draws in lgpfit (", N1, " vs. ", N2,
    ")! please report a bug!"
  )
  if (N1 != N2) {
    errors <- c(errors, msg)
  }

  # Validate postproc_results slot
  if (contains_postproc(object)) {
    ppn <- names(object@postproc_results)
    if (ppn != "pred") {
      errors <- c(errors, "allowed names for postproc_results: {pred}!")
    }
  }
  return_true_or_errors(errors)
}

return_true_or_errors <- function(errors) {
  if (length(errors) > 0) errors else TRUE
}

#' An S4 class to represent an lgp expression
#'
#' @slot covariate name of a covariate
#' @slot fun function name
#' @seealso See \code{\link{operations}} for performing arithmetics
#' on \linkS4class{lgprhs}, \linkS4class{lgpterm} and \linkS4class{lgpexpr}
#' objects.
lgpexpr <- setClass("lgpexpr",
  representation = representation(
    covariate = "character",
    fun = "character"
  ),
  prototype(covariate = "", fun = ""),
  validity = validate_lgpexpr
)

#' An S4 class to represent one formula term
#'
#' @slot factors a list of at most two \linkS4class{lgpexpr}s
#' @seealso See \code{\link{operations}} for performing arithmetics
#' on \linkS4class{lgprhs}, \linkS4class{lgpterm} and \linkS4class{lgpexpr}
#' objects.
lgpterm <- setClass("lgpterm", slots = c(factors = "list"))

#' An S4 class to represent the right-hand side of an lgp formula
#'
#' @slot summands a list of one or more \linkS4class{lgpterm}s
#' @seealso See \code{\link{operations}} for performing arithmetics
#' on \linkS4class{lgprhs}, \linkS4class{lgpterm} and \linkS4class{lgpexpr}
#' objects.
lgprhs <- setClass("lgprhs", slots = c(summands = "list"))

#' An S4 class to represent an lgp formula
#'
#' @slot terms an object of class \linkS4class{lgprhs}
#' @slot y_name name of the response variable
#' @slot call original formula call
#' @seealso See \code{\link{operations}} for performing arithmetics
#' on \linkS4class{lgprhs}, \linkS4class{lgpterm} and \linkS4class{lgpexpr}
#' objects.
lgpformula <- setClass("lgpformula",
  representation = representation(
    call = "character",
    y_name = "character",
    terms = "lgprhs"
  ),
  validity = validate_lgpformula
)

#' An S4 class to represent variable scaling
#'
#' @slot loc original location (mean)
#' @slot scale original scale (standard deviation)
#' @slot var_name variable name
lgpscaling <- setClass("lgpscaling",
  representation = representation(
    loc = "numeric",
    scale = "numeric",
    var_name = "character"
  ),
  prototype = prototype(loc = 1.0, scale = 1.0, var_name = "unknown"),
  validity = validate_lgpscaling
)


#' @rdname validate
validate_GaussianPrediction <- function(object) {
  errors <- c()
  f_comp_mean <- object@f_comp_mean
  f_comp_std <- object@f_comp_std
  errors <- c(errors, validate_lengths(f_comp_mean, f_comp_std))
  D1 <- dim(object@f_mean)
  D2 <- dim(object@f_std)
  D3 <- dim(object@y_mean)
  D4 <- dim(object@y_std)
  D <- list(D1, D2, D3, D4)
  L <- length(f_comp_mean)
  for (j in seq_len(L)) {
    m <- object@f_comp_mean[[j]]
    s <- object@f_comp_std[[j]]
    D <- c(D, list(dim(m), dim(s)))
  }
  errors <- c(errors, validate_dimension_list(D))
  if (is.null(object@x)) {
    errors <- c(errors, "x can't be NULL!")
  }
  out <- if (length(errors) > 0) errors else TRUE
  return(out)
}

#' @rdname validate
validate_Prediction <- function(object) {
  L <- length(object@f_comp)
  errors <- c()
  D1 <- dim(object@f)
  D2 <- dim(object@h)
  D <- list(D1, D2)
  for (j in seq_len(L)) {
    fj <- object@f_comp[[j]]
    D <- c(D, list(dim(fj)))
  }
  errors <- c(errors, validate_dimension_list(D))
  out <- if (length(errors) > 0) errors else TRUE
  return(out)
}

#' An S4 class to represent an additive GP model
#'
#' @slot formula An object of class \linkS4class{lgpformula}
#' @slot data The original unmodified data.
#' @slot stan_input The data to be given as input to \code{rstan::sampling}
#' @slot var_names List of variable names grouped by type.
#' @slot var_scalings A named list with fields
#' \itemize{
#'   \item \code{y} - Response variable normalization function and its
#'   inverse operation. Must be an \linkS4class{lgpscaling} object.
#'   \item \code{x_cont} - Continuous covariate normalization functions and
#'   their inverse operations. Must be a named list with each element is an
#'   \linkS4class{lgpscaling} object.
#' }
#' @slot var_info A named list with fields
#' \itemize{
#'   \item \code{x_cat_levels} - Names of the levels of categorical covariates
#'   before converting from factor to numeric.
#' }
#' @slot info Other info in text format.
#' @slot sample_f Whether the signal \code{f} is sampled or marginalized.
#' @slot full_prior Complete prior information.
#' @param object The object for which to apply a class method.
lgpmodel <- setClass("lgpmodel",
  representation = representation(
    model_formula = "lgpformula",
    data = "data.frame",
    stan_input = "list",
    var_names = "list",
    var_scalings = "list",
    var_info = "list",
    info = "list",
    sample_f = "logical",
    full_prior = "list"
  )
)

#' An S4 class to represent the output of the \code{lgp} function
#'
#' @slot stan_fit An object of class \code{stanfit}.
#' @slot model An object of class \code{lgpmodel}.
#' @slot num_draws Total number of parameter draws.
#' @slot postproc_results A named list containing possible postprocessing
#' results.
#' @param object The object for which to apply a class method.
#' @seealso For extracting parameter draws, see \code{\link{get_draws}},
#' or the \code{rstan} methods for \code{stanfit} objects.
lgpfit <- setClass("lgpfit",
  slots = c(
    stan_fit = "stanfit",
    model = "lgpmodel",
    num_draws = "numeric",
    postproc_results = "list"
  ),
  validity = validate_lgpfit
)

#' An S4 class to represent a data set simulated using the additive GP
#' formalism
#'
#' @slot data the actual data
#' @slot response name of the response variable in the data
#' @slot components the drawn function components
#' @slot kernel_matrices the covariance matrices for each gp
#' @slot info A list with fields
#' \itemize{
#'   \item \code{par_ell} the used lengthscale parameters
#'   \item \code{par_cont} the parameters used to generate the continuous
#'   covariates
#'   \item \code{p_signal} signal proportion
#' }
#' @slot effect_times A list with fields
#' \itemize{
#'   \item \code{true} possible true effect times that generate the disease
#'   effect
#'   \item \code{observed} possible observed effect times
#' }
lgpsim <- setClass("lgpsim",
  representation = representation(
    data = "data.frame",
    response = "character",
    components = "data.frame",
    kernel_matrices = "array",
    effect_times = "list",
    info = "list"
  )
)


#' An S4 class to represent analytically computed predictive distributions
#' (conditional on hyperparameters) of an additive GP model
#'
#' @slot f_comp_mean component means
#' @slot f_comp_std component standard deviations
#' @slot f_mean signal mean (on normalized scale)
#' @slot f_std signal standard deviation (on normalized scale)
#' @slot y_mean predictive mean (on original data scale)
#' @slot y_std predictive standard deviation (on original data scale)
#' @slot x a data frame of points (covariate values) where the
#' function posteriors or predictive distributions have been evaluated
#' @param object \linkS4class{GaussianPrediction} object for which to apply a
#' class method.
#' @seealso \linkS4class{Prediction}
GaussianPrediction <- setClass("GaussianPrediction",
  representation = representation(
    f_comp_mean = "list",
    f_comp_std = "list",
    f_mean = "matrix",
    f_std = "matrix",
    y_mean = "matrix",
    y_std = "matrix",
    x = "data.frame"
  ),
  validity = validate_GaussianPrediction
)

#' An S4 class to represent prior or posterior
#' draws from an additive function distribution.
#'
#' @slot f_comp component draws
#' @slot f signal draws
#' @slot h predictions (signal draws + scaling factor \code{c_hat},
#' transformed through inverse link function)
#' @slot x a data frame of points (covariate values) where the
#' functions/predictions have been evaluated/sampled
#' @slot extrapolated Boolean value telling if the function draws are
#' original MCMC draws or if they have been created by extrapolating
#' such draws.
#' @param object \linkS4class{Prediction} object for which to apply a class
#' method.
#' @seealso \linkS4class{GaussianPrediction}
Prediction <- setClass("Prediction",
  representation = representation(
    f_comp = "list",
    f = "matrix",
    h = "matrix",
    x = "data.frame",
    extrapolated = "logical"
  ),
  validity = validate_Prediction
)

#' An S4 class to represent input for kernel matrix computations
#'
#' @slot input Common input (for example parameter values).
#' @slot K_input Input for computing kernel matrices between data points
#' (\code{N} x \code{N}). A list.
#' @slot Ks_input Input for computing kernel matrices between data and output
#' points (\code{P} x \code{N}). A list.
#' @slot Kss_input Input for computing kernel matrices between output
#' points (\code{P} x \code{P}). A list, empty if \code{full_covariance=FALSE}.
#' @slot comp_names Component names (character vector).
#' @slot full_covariance Boolean value determining if this can compute
#' full predictive covariance matrices (or just marginal variance at each point).
#' @slot no_separate_output_points Boolean value determining if
#' \code{Ks_input} and \code{Kss_input} are the same thing. Using this
#' knowledge can reduce unnecessary computations of kernel matrices.
#' @slot STREAM external pointer (for calling 'Stan' functions)
#' @param object The object for which to call a class method.
KernelComputer <- setClass("KernelComputer",
  representation = representation(
    input = "list",
    K_input = "list",
    Ks_input = "list",
    Kss_input = "list",
    comp_names = "character",
    full_covariance = "logical",
    STREAM = "externalptr",
    no_separate_output_points = "logical"
  )
)

# Class info for show methods
class_info <- function(class_name) {
  str <- paste0(
    "An object of class ", class_name, ". See ?",
    class_name, " for more info."
  )
  return(str)
}

# Class info for show methods of function posterior objects
class_info_fp <- function(class_name, comp_names, D) {
  comp_str <- paste(comp_names, collapse = ", ")
  desc <- class_info(class_name)
  desc <- paste0(desc, "\n - ", D[1], " components: ", comp_str)
  desc <- paste0(desc, "\n - ", D[2], " parameter set(s)")
  desc <- paste0(desc, "\n - ", D[3], " evaluation points")
  return(desc)
}

# Check that all listed dimensions are equal
validate_dimension_list <- function(dims) {
  errors <- c()
  L <- length(dims)
  d1 <- dims[[1]]
  K <- length(d1)
  for (j in seq_len(L)) {
    dj <- dims[[j]]
    for (k in seq_len(K)) {
      if (dj[k] != d1[k]) {
        msg <- paste0(k, "th dimensions of elements 1 and ", j, " differ!")
        errors <- c(errors, msg)
      }
    }
  }
  return(errors)
}


# Check that arguments have equal lengths
validate_lengths <- function(a, b) {
  errors <- c()
  L1 <- length(a)
  L2 <- length(b)
  if (L1 != L2) {
    msg <- "lengths do not agree!"
    errors <- c(errors, msg)
  }
  return(errors)
}


# S4 GENERICS -------------------------------------------------------------

#' S4 generics for lgpfit, lgpmodel, and other objects
#'
#' @param object object for which to apply the generic
#' @param digits number of digits to show
#' @param ... additional optional arguments to pass
#' @name s4_generics
#' @seealso To find out which methods have been implemented for which classes,
#' see \linkS4class{lgpfit}, \linkS4class{lgpmodel},
#' \linkS4class{Prediction} and \linkS4class{GaussianPrediction}.
#' @return
#' \itemize{
#'    \item \code{parameter_info} returns a data frame with
#'    one row for each parameter and columns
#'    for parameter name, parameter bounds, and the assigned prior
#'    \item \code{component_info} returns a data frame with one row for
#'    each model component, and columns encoding information about
#'    model components
#'    \item \code{covariate_info} returns a list with names
#'    \code{continuous} and \code{categorical}, with information about
#'    both continuous and categorical covariates
#'    \item \code{component_names} returns a character vector with
#'    component names
#'    \item \code{get_model} for \linkS4class{lgpfit} objects
#'    returns an \linkS4class{lgpmodel}
#'    \item \code{is_f_sampled} returns a logical value
#'    \item \code{get_stanfit} returns a \code{stanfit} (rstan)
#'    \item \code{postproc} applies postprocessing and returns an
#'    updated \linkS4class{lgpfit}
#'    \item \code{clear_postproc} removes postprocessing information and
#'    returns an updated \linkS4class{lgpfit}
#'    \item \code{num_paramsets}, \code{num_evalpoints} and
#'    \code{num_components} return an integer
#'
#'
#' }
NULL

#' @describeIn s4_generics Get parameter information (priors etc.).
setGeneric(
  "parameter_info",
  function(object, digits) standardGeneric("parameter_info")
)

#' @describeIn s4_generics Get component information.
setGeneric(
  "component_info", function(object) standardGeneric("component_info")
)

#' @describeIn s4_generics Get covariate information.
setGeneric(
  "covariate_info", function(object) standardGeneric("covariate_info")
)

#' @describeIn s4_generics Get component names.
setGeneric(
  "component_names", function(object) standardGeneric("component_names")
)

#' @describeIn s4_generics Get \linkS4class{lgpmodel} object.
setGeneric(
  "get_model",
  function(object) standardGeneric("get_model")
)

#' @describeIn s4_generics Determine if signal f is sampled or marginalized.
setGeneric(
  "is_f_sampled",
  function(object) standardGeneric("is_f_sampled")
)

#' @describeIn s4_generics Extract stanfit object.
setGeneric(
  "get_stanfit",
  function(object) standardGeneric("get_stanfit")
)

#' @describeIn s4_generics Perform postprocessing.
setGeneric(
  "postproc", function(object, ...) standardGeneric("postproc")
)

#' @describeIn s4_generics Determine if object contains postprocessing
#' information.
setGeneric(
  "contains_postproc", function(object) standardGeneric("contains_postproc")
)

#' @describeIn s4_generics Clear postprocessing information (to reduce
#' size of object).
setGeneric(
  "clear_postproc", function(object) standardGeneric("clear_postproc")
)

#' @describeIn s4_generics Get number of parameter sets.
setGeneric(
  "num_paramsets", function(object) standardGeneric("num_paramsets")
)

#' @describeIn s4_generics Get number of points where posterior is evaluated.
setGeneric(
  "num_evalpoints", function(object) standardGeneric("num_evalpoints")
)

#' @describeIn s4_generics Get number of model components.
setGeneric(
  "num_components", function(object) standardGeneric("num_components")
)
