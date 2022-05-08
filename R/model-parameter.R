# Base class for parameters
Parameter <- R6::R6Class("Parameter",
  inherit = Printable,
  public = list(

    # constructor
    initialize = function(name, prior) {
      private$set_name(name)
      checkmate::assert_class(prior, "ContinuousDistribution")
      private$set_prior(prior)
    },

    # get prior
    get_prior = function() {
      private$prior
    }
  ),
  private = list(
    prior = NULL,

    # long notation
    notation_long = function() {
      paste0(
        private$name_string(), " ~ ", self$get_prior()$notation(FALSE)
      )
    },

    # short notation
    notation_short = function() {
      paste0(
        private$name_string(), " ~ ", self$get_prior()$notation(TRUE)
      )
    },

    # set prior
    set_prior = function(prior) {
      private$prior <- prior
    }
  )
)

# Base class for vector parameters
VectorParameter <- R6::R6Class("VectorParameter",
  inherit = Parameter,
  public = list(
    initialize = function(name, length) {
      private$set_name(name)
      private$set_length(length)
    },

    # get length
    get_length = function() {
      private$length
    }
  ),
  private = list(
    length = NULL,

    # set length
    set_length = function(length) {
      checkmate::assert_integerish(length, lower = 1)
      private$length <- length
    },

    # long notation
    notation_long = function() {
      paste0(
        "VectorParameter(", private$name_string(),
        ", length = ", number_string(self$get_length()), ")"
      )
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
