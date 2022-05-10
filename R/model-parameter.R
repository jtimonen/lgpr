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
    },

    # formatted name string
    name_string = function() {
      parameter_string(self$get_name())
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
        self$name_string(), "<lower = ", number_string(b$lower),
        ", upper = ", number_string(b$upper), ">",
        " ~ ", self$get_prior()$notation(FALSE)
      )
    },

    # short notation
    notation_short = function() {
      b <- self$get_bounds()
      paste0(
        self$name_string(), "<", number_string(b$lower), ", ",
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
          self$name_string(), "[", j, "]",
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
          self$name_string(), "[", j, "]",
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
