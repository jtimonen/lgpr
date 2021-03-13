return_true_or_errors <- function(errors) {
  if (length(errors) > 0) errors else TRUE
}

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
  msg <- paste0("invalid num_draws in lgpfit (", N1, " vs. ", N2 ,
                ")! please report a bug!")
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
#' @slot stan_dat The data to be given as input to \code{rstan::sampling}
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
#' @seealso
#' \code{\link{model_getters}}
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
#' @seealso For complete info on accessing the properties of the
#' \code{stan_fit} slot, see
#' \href{https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html}{here}.
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
#' @seealso
#' For visualizing, see \code{\link{plot_sim}}.
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
#' @seealso \linkS4class{FunctionPosterior}
FunctionDraws <- setClass("FunctionDraws",
  representation = representation(
    components = "data.frame",
    x = "data.frame",
    model = "lgpmodel"
  ),
  validity = validate_FunctionDraws
)
