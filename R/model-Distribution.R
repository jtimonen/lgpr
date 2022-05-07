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


# Log-Normal distribution
LogNormalDistribution <- R6::R6Class("LogNormalDistribution",
  inherit = ContinuousDistribution,
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
