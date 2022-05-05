# Base class for printable objects
Printable <- R6::R6Class("Printable",
  public = list(

    # constructor
    initialize = function(name) {
      private$set_name(name)
    },

    # print info
    print = function() {
      cat(private$notation_long())
    },

    # print more verbose info
    print_verbose = function() {
      cat(private$description())
    },

    # getter for name
    get_name = function() {
      private$name
    }
  ),
  private = list(
    name = NULL,

    # set and check name
    set_name = function(name) {
      checkmate::assert_string(name, min.chars = 1)
      private$name <- name
    },

    # description
    description = function() {
      classes <- setdiff(class(self), c("Printable", "R6"))
      str <- paste(classes, collapse = ", ")
      desc <- paste0("An object of class(es): ", str, "\n")
      desc <- paste0(desc, "Short notation: ", private$notation_short(), "\n")
      paste0(desc, "Long notation: ", private$notation_long(), "\n")
    },

    # formatted name string
    name_string = function() {
      variable_string(private$name)
    },

    # short notation
    notation_short = function() {
      private$name_string()
    },

    # long notation
    notation_long = function() {
      private$notation_short()
    }
  )
)


# Base class for variables
Variable <- R6::R6Class("Variable",
  inherit = Printable,
  public = list(),
  private = list(

    # long notation
    notation_long = function() {
      paste0("Variable(", private$name_string(), ")")
    }
  )
)

# Continuous variable
ContinuousVariable <- R6::R6Class("ContinuousVariable",
  inherit = Variable,
  public = list(
    initialize = function(name, mean = 0.0, sd = 1.0) {
      private$set_name(name)
      private$set_scale(mean, sd)
    },
    get_mean = function() {
      private$mean
    },
    get_sd = function() {
      private$sd
    }
  ),
  private = list(
    mean = NULL,
    sd = NULL,

    # set and check scaling mean and standard deviation
    set_scale = function(mean, sd) {
      checkmate::assert_numeric(mean)
      checkmate::assert_numeric(sd, lower = 1e-12)
      private$mean <- mean
      private$sd <- sd
    },

    # long notation
    notation_long = function() {
      paste0(
        "ContinuousVariable(", private$name_string(),
        ", mean = ", number_string(private$mean), ", sd = ",
        number_string(private$sd), ")"
      )
    }
  )
)

# Continuous variable
CategoricalVariable <- R6::R6Class("CategoricalVariable",
  inherit = Variable,
  public = list(
    initialize = function(name, categories = "foo") {
      private$set_name(name)
      private$set_categories(categories)
    },
    get_categories = function() {
      private$categories
    }
  ),
  private = list(
    categories = NULL,

    # set and check number of categories
    set_categories = function(categories) {
      num_unique_categories <- length(unique(categories))
      num_categories <- length(categories)
      checkmate::assert_true(num_unique_categories == num_categories)
      private$categories <- categories
    },

    # get number of categories
    num_categories = function() {
      length(private$categories)
    },

    # long notation
    notation_long = function() {
      paste0(
        "CategoricalVariable(", private$name_string(),
        ", num_categories = ", number_string(private$num_categories()), ")"
      )
    }
  )
)
