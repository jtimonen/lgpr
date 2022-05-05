#' Kernel description classes
#'
#' @section Methods:
#' \code{$desc()} Get kernel description.
#'
#' \code{$print()} Print the kernel description.
#'
#' \code{$parameters()} Get all \code{Parameter}s of the kernel as a list.
#'
#' \code{$variables()} Get all input \code{Variable}s of the kernel as a list.
#' @name Kernel-classes
NULL

Kernel <- R6::R6Class("Kernel",
  public = list(

    # constructor
    initialize = function() {
    },

    # description, must be overridden by inheriting class
    desc = function() {
      stop("Desc not implemented!")
    },

    # print info
    print = function() {
      cat(self$desc(), "\n")
    },

    # get all parameters as a list
    parameters = function() {
      list()
    },

    # get all input variables as a list
    variables = function() {
      list()
    }
  )
)

# Exponentiated quadratic kernel
ExpQuadKernel <- R6::R6Class("ExpQuadKernel",
  inherit = Kernel,
  public = list(
    x = NULL, # continuous input
    ell = NULL, # lengthscale parameter

    # constructor
    initialize = function(x, ell) {
      checkmate::assert_class(x, "ContinuousVariable")
      checkmate::assert_class(ell, "LengthscaleParameter")
      self$x <- x
      self$ell <- ell
    },

    # get all parameters as a list
    parameters = function() {
      list(self$ell)
    },

    # get all variables as a list
    variables = function() {
      list(self$x)
    },

    # description
    desc = function() {
      paste0("ExpQuadKernel(", self$x, ")")
    }
  )
)

# Input-warping kernel
InputWarpKernel <- R6::R6Class("InputWarpKernel",
  inherit = Kernel,
  public = list(
    variable = NULL,
    lengthscale = NULL,
    steepness = NULL,

    # constructor
    initialize = function(variable, lengthscale, steepness) {
      checkmate::assert_class(variable, "ContinuousVariable")
      checkmate::assert_class(lengthscale, "LengthscaleParameter")
      checkmate::assert_class(steepness, "WarpingSteepnessParameter")
      self$variable <- variable
      self$lengthscale <- lengthscale
      self$steepness <- steepness
    },

    # get all parameters as a list
    parameters = function() {
      list(self$lengthscale, self$steepness)
    },

    # description
    desc = function() {
      paste0("InputWarpKernel(", self$variable$name, ")")
    }
  )
)

# Categorical kernel
CategoricalKernel <- R6::R6Class("CategoricalKernel",
  inherit = Kernel,
  public = list(
    variable = NULL,

    # constructor
    initialize = function(variable) {
      checkmate::assert_class(variable, "CategoricalVariable")
      self$variable <- variable
    },

    # description
    desc = function() {
      paste0("CategoricalKernel(", self$variable$name, ")")
    }
  )
)

# Zero-sum kernel
ZeroSumKernel <- R6::R6Class("ZeroSumKernel",
  inherit = Kernel,
  public = list(
    variable = NULL,

    # constructor
    initialize = function(variable) {
      checkmate::assert_class(variable, "CategoricalVariable")
      self$variable <- variable
    },

    # description
    desc = function() {
      paste0("ZeroSumKernel(", self$variable$name, ")")
    }
  )
)


# Heterogeneity kernel
HeterogeneityKernel <- R6::R6Class("HeterogeneityKernel",
  inherit = Kernel,
  public = list(
    variable = NULL,
    beta = NULL,

    # constructor
    initialize = function(variable, beta) {
      checkmate::assert_class(variable, "CategoricalVariable")
      checkmate::assert_class(beta, "ParameterVector")
      checkmate::assert_class(beta$param, "UnitIntervalParameter")
      checkmate::assert_true(beta$length == variable$num_categories())
      self$variable <- variable
      self$beta <- beta
    },

    # get all parameters as a list
    parameters = function() {
      list(self$params)
    },

    # description
    desc = function() {
      paste0("HeterogeneityKernel(", self$variable$name, ")")
    }
  )
)
