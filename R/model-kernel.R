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
  inherit = Printable,
  public = list(

    # get all parameters as a list
    get_parameters = function() {
      list()
    },

    # get all input variables as a list
    get_variables = function() {
      list()
    }
  ),
  private = list(
    # short notation
    notation_short = function() {
      s <- lapply(self$get_variables(), function(x) x$name_string())
      str1 <- paste(s, collapse = ", ")
      s <- lapply(self$get_parameters(), function(x) x$name_string())
      str2 <- paste(s, collapse = ", ")
      paste0(
        self$name_string(), "(", str1, " | ", str2, ")"
      )
    },

    # long notation
    notation_long = function() {
      s <- lapply(self$get_variables(), function(x) x$notation(TRUE))
      str1 <- paste(s, collapse = ", ")
      s <- lapply(self$get_parameters(), function(x) x$notation(TRUE))
      str2 <- paste(s, collapse = "\n  ")
      paste0(
        "Kernel name: ", self$name_string(), "\nVariables: ", str1,
        "\nParameters: {\n  ", str2, "\n}"
      )
    }
  )
)

# Exponentiated quadratic kernel
ExpQuadKernel <- R6::R6Class("ExpQuadKernel",
  inherit = Kernel,
  public = list(

    # constructor
    initialize = function(name, x, ell_prior) {
      private$set_name(paste0("k_", name))
      checkmate::assert_class(x, "ContinuousVariable")
      checkmate::assert_class(ell_prior, "ContinuousDistribution")
      private$x <- x
      ell_name <- paste0("ell_", name)
      private$ell <- Parameter$new(name = ell_name, prior = ell_prior)
    },

    # get all parameters as a list
    get_parameters = function() {
      list(ell = private$ell)
    },

    # get all variables as a list
    get_variables = function() {
      list(x = private$x)
    }
  ),
  private = list(
    x = NULL, # continuous input
    ell = NULL # lengthscale parameter
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
