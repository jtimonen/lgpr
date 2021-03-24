# MAIN DOCUMENTATION PAGE -------------------------------------------------

#' The 'lgpr' package.
#'
#' @description Interpretable nonparametric modeling of longitudinal data
#' using additive Gaussian process regression. Contains functionality
#' for inferring covariate effects and assessing covariate relevances.
#' Models are specified using a convenient formula syntax, and can include
#' shared, group-specific, non-stationary, heterogeneous and temporally
#' uncertain effects. Bayesian inference for model parameters is performed
#' using Stan (\code{\link[rstan]{rstan}}).
#'
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords Gaussian processes, longitudinal data, Stan, covariate relevances,
#' interpretable models
#'
#'
#' @section Specifying and fitting models:
#' Core functionality of the package consists of creating and fitting an
#' additive GP model. Main functions are:
#' \itemize{
#'  \item \code{\link{lgp}}: Specify and fit a model with one command.
#'  \item \code{\link{create_model}}: Define an \linkS4class{lgpmodel} object.
#'  \item \code{\link{sample_model}}: Sample model parameters and create an
#'  an \linkS4class{lgpfit} object.
#'  \item \code{\link{optimize_model}} (experimental): Optimize model
#'  parameters.
#' }
#'
#' @section Data:
#' The data that you wish to analyze with \code{lgpr} should be in an R
#' \code{data.frame} where columns correspond to measured variables and rows
#' corresponds to obervations. Some functions that can help working with such
#' data frames are:
#' \itemize{
#'  \item \code{\link{plot_data}}: Visualizing data.
#'  \item \code{\link{new_x}}: Creating new test points where the posterior
#'  distribution of any function component or sum of all components, or the
#'  posterior predictive distribution can be computed after model fitting.
#'  \item Other functions: \code{\link{add_factor}},
#'  \code{\link{add_factor_crossing}}, \code{\link{add_dis_age}},
#'  \code{\link{adjusted_c_hat}}.
#' }
#'
#' @section Tutorials and case studies:
#' See \url{https://jtimonen.github.io/lgpr-usage/index.html}.
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
#'   \href{https://doi.org/10.1093/bioinformatics/btab021}{url}.
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
  N1 <- nrow(get_draws(object, NULL, NULL, pars = c("alpha")))
  N2 <- object@num_draws
  msg <- paste0(
    "invalid num_draws in lgpfit (", N1, " vs. ", N2,
    ")! please report a bug!"
  )
  if (N1 != N2) {
    errors <- c(errors, msg)
  }
  return_true_or_errors(errors)
}

#' @rdname validate
validate_FunctionPosterior <- function(object) {
  errors <- c()
  cn1 <- colnames(object@components)
  cn2 <- colnames(object@total)
  tgt1 <- c("paramset", "component", "eval_point", "mean", "std")
  tgt2 <- c("paramset", "eval_point", "mean", "std", "sigma")

  if (!all.equal(cn1, tgt1)) {
    errors <- c(errors, "invalid components data frame!")
  }
  if (!all.equal(cn2, tgt2)) {
    errors <- c(errors, "invalid total data frame!")
  }
  return_true_or_errors(errors)
}

#' @rdname validate
validate_FunctionDraws <- function(object) {
  errors <- c()
  cn <- colnames(object@components)
  tgt <- c("paramset", "component", "eval_point", "value")
  if (!all.equal(cn, tgt)) {
    errors <- c(errors, "invalid components data frame!")
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

#' An S4 class to represent an lgp model
#'
#' @slot formula An object of class \linkS4class{lgpformula}
#' @slot data The original unmodified data.
#' @slot stan_input The data to be given as input to \code{rstan::sampling}
#' @slot var_names List of variable names grouped by type.
#' @slot var_scalings A named list with fields
#' \itemize{
#'   \item \code{y} - Response variable normalization function and its
#'   inverse operation. Must be an \linkS4class{lgpscaling} object.
#'   \item \code{x_cont} - Continous covariate normalization functions and
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
#' @param object \linkS4class{lgpmodel} object for which to apply a class
#' method.
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
#' @param object \linkS4class{lgpfit} object for which to apply a class method.
#' @param ... optional arguments passed to a subroutine
lgpfit <- setClass("lgpfit",
  slots = c(
    stan_fit = "stanfit",
    model = "lgpmodel",
    num_draws = "numeric"
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

#' An S4 class to represent posterior distributions of an additive marginal GP
#'
#' @slot components Data frame representing the posterior distribution of
#' each model component (on normalized scale).
#' @slot total Data frame representing the posterior distribution of the sum of
#' the components (on normalized scale).
#' @slot x The evaluation points (values of covariates) where the posteriors
#' are evaluated.
#' @slot model The \linkS4class{lgpmodel} for which these posteriors are
#' computed. Contains important information about how to scale the total
#' posterior from normalized scale to the original scale.
#' @param object \linkS4class{FunctionPosterior} object for which to apply a
#' class method.
#' @seealso \linkS4class{FunctionDraws}
FunctionPosterior <- setClass("FunctionPosterior",
  representation = representation(
    components = "data.frame",
    total = "data.frame",
    x = "data.frame",
    model = "lgpmodel"
  ),
  validity = validate_FunctionPosterior
)

#' An S4 class to represent posterior or prior function (component)
#' draws of an additive latent GP model
#'
#' @slot components Data frame representing the draws of each model component.
#' @slot x The evaluation points (values of covariates) where the function
#' draws are evaluated.
#' @slot model The \linkS4class{lgpmodel} for which these draws are
#' computed. Contains important information about how to transform the
#' total draws through an inverse link function.
#' @param object \linkS4class{FunctionDraws} object for which to apply a class
#' method.
#' @seealso \linkS4class{FunctionPosterior}
FunctionDraws <- setClass("FunctionDraws",
  representation = representation(
    components = "data.frame",
    x = "data.frame",
    model = "lgpmodel"
  ),
  validity = validate_FunctionDraws
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
class_info_fp <- function(class_name, num_comps, D) {
  desc <- class_info(class_name)
  desc <- paste0(desc, "\n - ", num_comps, " components")
  desc <- paste0(desc, "\n - ", D[1], " parameter sets")
  desc <- paste0(desc, "\n - ", D[2], " evaluation points")
  return(desc)
}


# S4 GENERICS -------------------------------------------------------------

setGeneric(
  "parameter_info",
  function(object, digits) standardGeneric("parameter_info")
)

setGeneric(
  "component_info", function(object) standardGeneric("component_info")
)

setGeneric(
  "covariate_info", function(object) standardGeneric("covariate_info")
)

setGeneric(
  "component_names", function(object) standardGeneric("component_names")
)

setGeneric(
  "get_model",
  function(object) standardGeneric("get_model")
)

setGeneric(
  "get_stanfit",
  function(object) standardGeneric("get_stanfit")
)

setGeneric(
  "get_draws",
  function(object, draws = NULL, reduce = NULL, ...) standardGeneric("get_draws")
)
