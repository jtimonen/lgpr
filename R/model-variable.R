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
    },

    # get notation
    notation = function(short = FALSE) {
      if (short) {
        return(private$notation_short())
      }
      private$notation_long()
    },

    # formatted name string
    name_string = function() {
      self$get_name()
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
      desc <- paste0(desc, "Short notation: ", self$notation(TRUE), "\n")
      paste0(desc, "Long notation: ", self$notation(FALSE), "\n")
    },

    # short notation
    notation_short = function() {
      self$name_string()
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
  public = list(
    # formatted name string
    name_string = function() {
      variable_string(self$get_name())
    }
  ),
  private = list(

    # long notation
    notation_long = function() {
      paste0("Variable(", self$name_string(), ")")
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
      checkmate::assert_numeric(sd, lower = small_number())
      private$mean <- mean
      private$sd <- sd
    },

    # long notation
    notation_long = function() {
      paste0(
        "ContinuousVariable(", self$name_string(),
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
        "CategoricalVariable(", self$name_string(),
        ", num_categories = ", number_string(private$num_categories()), ")"
      )
    }
  )
)
