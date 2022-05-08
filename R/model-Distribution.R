# Base class for distribution
Distribution <- R6::R6Class("Distribution",
  inherit = Printable,
  public = list(
    initialize = function() {
    },

    # returns a named list of hyperparameter values
    get_hyperparams = function() {
      list()
    },
    get_name = function() {
      "Distribution"
    }
  ),
  private = list(

    # Get string that tells hyperparam names and values
    hyperparams_string = function(show_names = TRUE) {
      params <- self$get_hyperparams()
      str <- ""
      L <- length(params)
      plist <- list()
      for (j in seq_len(L)) {
        s <- ""
        if (show_names) {
          s <- paste0(names(params)[j], " = ")
        }
        plist[[j]] <- paste0(s, number_string(params[[j]]))
      }
      paste(plist, collapse = ", ")
    }
  )
)

# Base class for continuous distributions
ContinuousDistribution <- R6::R6Class("ContinuousDistribution",
  inherit = Distribution
)

# Base class for positive continuous distributions
PositiveContinuousDistribution <- R6::R6Class("PositiveContinuousDistribution",
  inherit = ContinuousDistribution
)

# Base class for unit interval distributions
UnitIntervalDistribution <- R6::R6Class("UnitIntervalDistribution",
  inherit = ContinuousDistribution
)

# Base class for discrete distributions
DiscreteDistribution <- R6::R6Class("DiscreteDistribution",
  inherit = Distribution
)

# Uniform distribution
Uniform <- R6::R6Class("Uniform",
  inherit = ContinuousDistribution,
  private = list(

    # short notation
    notation_short = function() {
      paste0("U")
    },

    # long notation
    notation_long = function() {
      paste0("Uniform")
    }
  )
)

# Normal distribution
NormalDistribution <- R6::R6Class("NormalDistribution",
  inherit = ContinuousDistribution,
  public = list(
    initialize = function(mean, sd) {
      checkmate::assert_numeric(mean, max.len = 1)
      checkmate::assert_numeric(sd, max.len = 1, lower = small_number())
      private$mean <- mean
      private$sd <- sd
    },
    get_hyperparams = function() {
      list(
        mean = private$mean,
        sd = private$sd
      )
    }
  ),
  private = list(
    mean = NULL,
    sd = NULL,

    # short notation
    notation_short = function() {
      paste0("N(", private$hyperparams_string(show_names = FALSE), ")")
    },

    # long notation
    notation_long = function() {
      paste0("Normal(", private$hyperparams_string(), ")")
    }
  )
)


# Student's t-distribution
StudentTDistribution <- R6::R6Class("StudentTDistribution",
  inherit = ContinuousDistribution,
  public = list(
    initialize = function(nu) {
      checkmate::assert_numeric(nu, max.len = 1, lower = small_number())
      private$nu <- nu
    },
    get_hyperparams = function() {
      list(
        nu = private$nu
      )
    }
  ),
  private = list(
    nu = NULL,

    # short notation
    notation_short = function() {
      paste0("t(", private$hyperparams_string(show_names = FALSE), ")")
    },

    # long notation
    notation_long = function() {
      paste0("Student-t(", private$hyperparams_string(), ")")
    }
  )
)

# Log-Normal distribution
LogNormalDistribution <- R6::R6Class("LogNormalDistribution",
  inherit = PositiveContinuousDistribution,
  public = list(
    initialize = function(logmean, logsd) {
      checkmate::assert_numeric(logmean, max.len = 1)
      checkmate::assert_numeric(logsd, max.len = 1, lower = small_number())
      private$logmean <- logmean
      private$logsd <- logsd
    },
    get_hyperparams = function() {
      list(
        logmean = private$logmean,
        logsd = private$logsd
      )
    }
  ),
  private = list(
    logmean = NULL,
    logsd = NULL,

    # short notation
    notation_short = function() {
      paste0("Log-N(", private$hyperparams_string(show_names = FALSE), ")")
    },

    # long notation
    notation_long = function() {
      paste0("LogNormal(", private$hyperparams_string(), ")")
    }
  )
)

# Gamma distribution
GammaDistribution <- R6::R6Class("GammaDistribution",
  inherit = PositiveContinuousDistribution,
  public = list(
    initialize = function(shape, inv_scale) {
      checkmate::assert_numeric(shape, max.len = 1, lower = small_number())
      checkmate::assert_numeric(inv_scale, max.len = 1, lower = small_number())
      private$shape <- shape
      private$inv_scale <- inv_scale
    },
    get_hyperparams = function() {
      list(
        shape = private$shape,
        inv_scale = private$inv_scale
      )
    }
  ),
  private = list(
    shape = NULL,
    inv_scale = NULL,

    # short notation
    notation_short = function() {
      paste0("Gamma(", private$hyperparams_string(show_names = FALSE), ")")
    },

    # long notation
    notation_long = function() {
      paste0("Gamma(", private$hyperparams_string(), ")")
    }
  )
)

# Inverse Gamma distribution
InverseGammaDistribution <- R6::R6Class("InverseGammaDistribution",
  inherit = PositiveContinuousDistribution,
  public = list(
    initialize = function(shape, scale) {
      checkmate::assert_numeric(shape, max.len = 1, lower = small_number())
      checkmate::assert_numeric(scale, max.len = 1, lower = small_number())
      private$shape <- shape
      private$scale <- scale
    },
    get_hyperparams = function() {
      list(
        shape = private$shape,
        scale = private$scale
      )
    }
  ),
  private = list(
    shape = NULL,
    scale = NULL,

    # short notation
    notation_short = function() {
      paste0("InvGamma(", private$hyperparams_string(show_names = FALSE), ")")
    },

    # long notation
    notation_long = function() {
      paste0("InverseGamma(", private$hyperparams_string(), ")")
    }
  )
)


# Beta distribution
BetaDistribution <- R6::R6Class("BetaDistribution",
  inherit = UnitIntervalDistribution,
  public = list(
    initialize = function(shape1, shape2) {
      checkmate::assert_numeric(shape1, max.len = 1, lower = small_number())
      checkmate::assert_numeric(shape2, max.len = 1, lower = small_number())
      private$shape1 <- shape1
      private$shape2 <- shape2
    },
    get_hyperparams = function() {
      list(
        shape1 = private$shape1,
        shape2 = private$shape2
      )
    }
  ),
  private = list(
    shape1 = NULL,
    shape2 = NULL,

    # short notation
    notation_short = function() {
      paste0("Beta(", private$hyperparams_string(show_names = FALSE), ")")
    },

    # long notation
    notation_long = function() {
      paste0("Beta(", private$hyperparams_string(), ")")
    }
  )
)
