# Base class for variables
Variable <- R6::R6Class("Variable",
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
      paste0("Variable(", self$name, ")")
    },

    # print info
    print = function() {
      cat(self$desc(), "\n")
    }
  )
)

# Continuous variable
ContinuousVariable <- R6::R6Class("ContinuousVariable",
  inherit = Variable,
  public = list(
    mean = NULL,
    sd = NULL,
    initialize = function(name, mean = 0.0, sd = 1.0) {
      self$set_name(name)
      self$set_scale(mean, sd)
    },

    # set and check scaling mean and standard deviation
    set_scale = function(mean, sd) {
      checkmate::assert_numeric(mean)
      checkmate::assert_numeric(sd, lower = 1e-12)
      self$mean <- mean
      self$sd <- sd
    }
  )
)

# Continuous variable
CategoricalVariable <- R6::R6Class("CategoricalVariable",
  inherit = Variable,
  public = list(
    categories = NULL,
    initialize = function(name, categories = "foo") {
      self$set_name(name)
      self$set_categories(categories)
    },

    # set and check number of categories
    set_categories = function(categories) {
      checkmate::assert_set_equal(unique(categories), categories)
      self$categories <- categories
    },

    # get number of categories
    num_categories = function() {
      length(self$categories)
    }
  )
)
