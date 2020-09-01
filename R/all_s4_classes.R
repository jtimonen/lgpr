#' @include s4-validate.R
NULL

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
  validity = check_lgpexpr
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
  validity = check_lgpformula
)

#' An S4 class to represent variable scaling and its inverse
#'
#' @slot fun function that performs normalization
#' @slot fun_inv inverse function of \code{fun}
#' @slot var_name variable name
lgpscaling <- setClass("lgpscaling",
  representation = representation(
    fun = "function",
    fun_inv = "function",
    var_name = "character"
  ),
  prototype = prototype(
    fun = function(x) {
      x
    },
    fun_inv = function(x) {
      x
    },
    var_name = "unknown"
  ),
  validity = check_lgpscaling
)

#' An S4 class to represent an lgp model
#'
#' @slot formula An object of class \linkS4class{lgpformula}
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
#' @slot stan_model_name Name of Stan model.
#' @section Related methods and functions:
#' \itemize{
#'     \item \code{\link{model_summary}}
#'     \item \code{\link{prior_summary}} and \code{\link{param_summary}}
#'     \item \code{\link{model_getters}}
#' }
lgpmodel <- setClass("lgpmodel",
  representation = representation(
    model_formula = "lgpformula",
    stan_input = "list",
    var_names = "list",
    var_scalings = "list",
    var_info = "list",
    info = "list",
    stan_model_name = "character"
  )
)

#' An S4 class to represent the output of the \code{lgp} function
#'
#' @slot stan_fit An object of class \code{stanfit}.
#' @slot model An object of class \code{lgpmodel}.
#' @section Related methods and functions:
#' \itemize{
#'     \item \code{\link{fit_summary}}
#'     \item \code{\link{get_draws}}
#'     \item \code{\link{get_posterior_f}}
#'     \item \code{\link{plot_posterior}}
#'     \item \code{\link{plot_posterior_warp}}
#'     \item In addition, all methods that work on \linkS4class{lgpmodel}
#'     objects work also on \linkS4class{lgpfit} objects.
#' }
#' @seealso For complete info on accessing the properties of the \code{stan_fit}
#' slot, see \href{https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html}{here}.
lgpfit <- setClass("lgpfit",
  slots = c(
    stan_fit = "stanfit",
    model = "lgpmodel"
  )
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
#' @section Related methods and functions:
#' \itemize{
#'     \item \code{\link{sim_plot}}
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
