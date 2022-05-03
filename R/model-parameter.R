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
CovariateUncertaintyParameter <- R6::R6Class("CovariateUncertaintyParameter",
  inherit = IntervalParameter,
  public = list()
)

# Effect heterogeneity parameter
EffectHeterogeneityParameter <- R6::R6Class("EffectHeterogeneityParameter",
  inherit = IntervalParameter,
  public = list(
    # constructor
    initialize = function(name) {
      self$set_name(name)
      self$lower <- 0
      self$upper <- 1
    }
  )
)
