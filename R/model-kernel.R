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
    }
  )
)

# Exponentiated quadratic kernel
ExpQuadKernel <- R6::R6Class("ExpQuadKernel",
  inherit = Kernel,
  public = list(
    variable = NULL,

    # constructor
    initialize = function(variable) {
      checkmate::assert_class(variable, "ContinuousVariable")
      self$variable <- variable
    },

    # description
    desc = function() {
      paste0("ExpQuadKernel(", self$variable$name, ")")
    }
  )
)
