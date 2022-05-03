# Base class for parameters
Parameter <- R6::R6Class("Parameter",
  public = list(
    name = NULL,

    # constructor
    initialize = function(name) {
      self$set_name(name)
    },

    # set and check name
    set_name = function(name) {
      checkmate::assert_string(name, min.chars = 1)
      self$name <- name
    },

    # description
    desc = function() {
      paste0("Parameter(", self$name, ")")
    },

    # print info
    print = function() {
      cat(self$desc(), "\n")
    }
  )
)

# Parameter that is restricted to be positive
PositiveParameter <- R6::R6Class("PositiveParameter",
  inherit = Parameter,
  public = list()
)

# Parameter that is restricted to be on some interval
IntervalParameter <- R6::R6Class("IntervalParameter",
  inherit = Parameter,
  public = list(
    upper = NULL,
    lower = NULL,

    # constructor
    initialize = function(name, lower = 0, upper = 1) {
      self$set_name(name)
      checkmate::assert_numeric(lower)
      checkmate::assert_numeric(upper)
      self$lower <- lower
      self$upper <- upper
    }
  )
)

# Parameter that is restricted to be on unit interval
UnitIntervalParameter <- R6::R6Class("UnitIntervalParameter",
  inherit = Parameter,
  public = list()
)

# Lengthscale parameter
LengthscaleParameter <- R6::R6Class("LengthscaleParameter",
  inherit = PositiveParameter,
  public = list()
)

# Warping steepness parameter
WarpingSteepnessParameter <- R6::R6Class("WarpingSteepnessParameter",
  inherit = PositiveParameter,
  public = list()
)

# Lengthscale parameter
MagnitudeParameter <- R6::R6Class("MagnitudeParameter",
  inherit = PositiveParameter,
  public = list()
)

# Covariate uncertainty parameter
ParameterVector <- R6::R6Class("ParameterVector",
  public = list(
    param = NULL,
    length = NULL,
    # constructor
    initialize = function(param, length) {
      checkmate::assert_class(param, "Parameter")
      checkmate::assert_integerish(length, lower = 1)
      self$param <- param
      self$length <- length
    },

    # description
    desc = function() {
      paste0(
        "ParameterVector(", self$param$desc(), ", ", self$length, ")"
      )
    },

    # print info
    print = function() {
      cat(self$desc(), "\n")
    }
  )
)
