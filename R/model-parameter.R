# Base class for parameters
Parameter <- R6::R6Class("Parameter",
  inherit = Printable,
  public = list(

    # constructor
    initialize = function(name, prior, lower = -Inf, upper = Inf) {
      private$set_name(name)
      private$set_prior(prior)
      private$set_bounds(lower, upper)
    },

    # get prior
    get_prior = function() {
      private$prior
    },

    # get bounds
    get_bounds = function() {
      list(
        lower = private$lower,
        upper = private$upper
      )
    }
  ),
  private = list(
    prior = NULL,
    lower = NULL,
    upper = NULL,

    # long notation
    notation_long = function() {
      b <- self$get_bounds()
      paste0(
        private$name_string(), "<lower = ", number_string(b$lower),
        ", upper = ", number_string(b$upper), ">",
        " ~ ", self$get_prior()$notation(FALSE)
      )
    },

    # short notation
    notation_short = function() {
      b <- self$get_bounds()
      paste0(
        private$name_string(), "<", number_string(b$lower), ", ",
        number_string(b$upper), ">",
        " ~ ", self$get_prior()$notation(TRUE)
      )
    },

    # set prior
    set_prior = function(prior) {
      checkmate::assert_class(prior, "ContinuousDistribution")
      private$prior <- prior
    },

    # set bounds
    set_bounds = function(lower, upper) {
      checkmate::assert_number(lower, finite = FALSE)
      checkmate::assert_number(upper, finite = FALSE)
      checkmate::assert_true(lower < upper)
      private$lower <- lower
      private$upper <- upper
    }
  )
)

# Base class for vector parameters
VectorParameter <- R6::R6Class("VectorParameter",
  inherit = Parameter,
  public = list(

    # constructor
    initialize = function(name, length, prior, lower = -Inf, upper = Inf) {
      checkmate::assert_integerish(length, lower = 1)
      private$length <- length
      private$set_name(name)
      private$set_prior(prior)
      private$set_bounds(lower, upper)
    },

    # get vector length
    get_length = function() {
      private$length
    },

    # get prior of parameter with index idx
    get_prior = function(idx) {
      L <- self$get_length()
      checkmate::assert_integerish(idx, lower = 1, upper = L)
      prior <- private$prior
      if (is(prior, "ContinuousDistribution")) {
        return(prior)
      }
      prior[[idx]]
    },

    # get bounds of parameter with index idx
    get_bounds = function(idx) {
      L <- self$get_length()
      checkmate::assert_integerish(idx, lower = 1, upper = L)
      lower <- private$lower
      upper <- private$upper
      if (length(lower) > 1) {
        lower <- lower[[idx]]
        upper <- upper[[idx]]
      }
      list(
        lower = lower,
        upper = upper
      )
    }
  ),
  private = list(
    length = NULL,
    prior = NULL,
    lower = NULL,
    upper = NULL,

    # long notation
    notation_long = function() {
      desc <- ""
      L <- self$get_length()
      for (j in seq_len(L)) {
        b <- self$get_bounds(idx = j)
        pr <- self$get_prior(idx = j)
        desc <- paste0(
          desc, "\n",
          private$name_string(), "[", j, "]",
          "<lower = ", number_string(b$lower), ", upper = ",
          number_string(b$upper), ">",
          " ~ ", pr$notation(FALSE)
        )
      }
      return(desc)
    },

    # short notation
    notation_short = function() {
      desc <- ""
      L <- self$get_length()
      for (j in seq_len(L)) {
        b <- self$get_bounds(idx = j)
        pr <- self$get_prior(idx = j)
        desc <- paste0(
          desc, "\n",
          private$name_string(), "[", j, "]",
          "<", number_string(b$lower), ", ",
          number_string(b$upper), ">",
          " ~ ", pr$notation(TRUE)
        )
      }
      return(desc)
    },

    # set prior
    set_prior = function(prior) {
      if (is(prior, "ContinuousDistribution")) {
        # do nothing
      } else {
        L <- self$get_length()
        checkmate::assert_list(prior, len = L)
        for (j in seq_len(L)) {
          checkmate::assert_class(prior[[j]], "ContinuousDistribution")
        }
      }
      private$prior <- prior
    },

    # set bounds
    set_bounds = function(lower, upper) {
      if (length(lower) == 1) {
        checkmate::assert_number(lower, finite = FALSE)
        checkmate::assert_number(upper, finite = FALSE)
        checkmate::assert_true(lower < upper)
      } else {
        L <- self$get_length()
        checkmate::assert_numeric(lower, finite = FALSE, len = L)
        checkmate::assert_numeric(upper, finite = FALSE, len = L)
        checkmate::assert_true(all(lower < upper))
      }
      private$lower <- lower
      private$upper <- upper
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
