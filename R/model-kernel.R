# Base class for kernels
Kernel <- R6::R6Class("Kernel",
  public = list(

    # constructor
    initialize = function() {
    },

    # description
    desc = function() {
      paste0("Kernel")
    },

    # print info
    print = function() {
      cat(self$desc(), "\n")
    },

    # get all parameters as a list
    parameters = function() {
      list()
    }
  )
)

# Exponentiated quadratic kernel
ExpQuadKernel <- R6::R6Class("ExpQuadKernel",
  inherit = Kernel,
  public = list(
    variable = NULL,
    lengthscale = NULL,

    # constructor
    initialize = function(variable, lengthscale) {
      checkmate::assert_class(variable, "ContinuousVariable")
      checkmate::assert_class(lengthscale, "LengthscaleParameter")
      self$variable <- variable
      self$lengthscale <- lengthscale
    },

    # get all parameters as a list
    parameters = function() {
      list(self$lengthscale)
    },

    # description
    desc = function() {
      paste0("ExpQuadKernel(", self$variable$name, ")")
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
