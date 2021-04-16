# Argument validation
#
# @param arg An argument to check.
# @return \code{TRUE} if the check passes, throw error otherwise.

# checks that argument has expected type
check_type <- function(arg, type) {
  arg_name <- deparse(substitute(arg))
  is_type <- is(arg, type)
  if (!is_type) {
    found <- as.character(class(arg))
    stop(
      "'", arg_name, "' must be an object of type <", type, ">!",
      " Found = <", found, ">",
      sep = ""
    )
  }
  TRUE
}

# checks that argument is numeric
check_numeric <- function(arg) {
  if (!is.numeric(arg)) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must be numeric! found = ", arg)
    stop(msg)
  }
  TRUE
}

# checks that argument is positive
check_positive <- function(arg) {
  check_numeric(arg)
  if (arg <= 0) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must be positive! found = ", arg)
    stop(msg)
  }
  TRUE
}

# checks that argument is non-negative
check_non_negative <- function(arg) {
  check_numeric(arg)
  if (arg < 0) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must be non-negative! found = ", arg)
    stop(msg)
  }
  TRUE
}

# checks that argument has only positive values
check_positive_all <- function(arg) {
  check_numeric(arg)
  if (any(arg <= 0)) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must have only positive values")
    stop(msg)
  }
  TRUE
}

# checks that argument has only non-negative values
check_non_negative_all <- function(arg) {
  check_numeric(arg)
  if (any(arg < 0)) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must have only non-negative values")
    stop(msg)
  }
  TRUE
}

# checks that argument has only integer values
check_integer_all <- function(arg) {
  check_numeric(arg)
  if (sum(round(arg) != arg) > 0) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must have only integer values")
    stop(msg)
  }
  TRUE
}

# checks checks that value is in some interval (strictly)
check_interval <- function(arg, lower, upper) {
  check_numeric(arg)
  if (arg < lower || arg > upper) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0(
      "<", arg_name, "> must be on the interval [",
      lower, ", ", upper, "]! found = ", arg
    )
    stop(msg)
  }
  TRUE
}

# checks that argument is not NULL
check_not_null <- function(arg) {
  arg_name <- deparse(substitute(arg))
  if (is.null(arg)) {
    stop(arg_name, " should not be NULL!")
  }
  TRUE
}

# checks checks that argument is NULL
# @param reason explanation why the argument should be NULL
check_null <- function(arg, reason = NULL) {
  arg_name <- deparse(substitute(arg))
  if (!is.null(arg)) {
    msg <- paste0("<", arg_name, "> should be NULL!")
    if (!is.null(reason)) {
      msg <- paste0(msg, " Reason: ", reason)
    }
    stop(msg)
  }
  TRUE
}

# checks that argument has expected length
check_length <- function(arg, expected) {
  L <- length(arg)
  ok <- (L == expected)
  if (!ok) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0(
      arg_name, " has length ", L, ", but its length should be ",
      expected, "!"
    )
    stop(msg)
  }
  TRUE
}

# checks that argument has length 1 or required other value
check_length_1_or <- function(arg, expected) {
  L <- length(arg)
  ok <- (L == expected) || (L == 1)
  if (!ok) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0(
      arg_name, " has length ", L, ", but its length should be ",
      expected, " or one!"
    )
    stop(msg)
  }
  TRUE
}

# checks that argument has at least min_length length
check_length_geq <- function(arg, min_length) {
  L <- length(arg)
  bad <- (L < min_length)
  if (bad) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0(
      arg_name, " has length ", L, ", but its length should be at least ",
      min_length, "!"
    )
    stop(msg)
  }
  TRUE
}

# checks that arguments a and b have same length
check_lengths <- function(a, b) {
  L1 <- length(a)
  L2 <- length(b)
  a_name <- deparse(substitute(a))
  b_name <- deparse(substitute(b))
  if (L1 != L2) {
    msg <- paste0(
      "lengths of ", a_name, " and ", b_name, " must match!",
      " found = ", L1, " and ", L2
    )
    stop(msg)
  }
  TRUE
}

# checks that argument is FALSE
check_false <- function(arg) {
  arg_name <- deparse(substitute(arg))
  check_not_null(arg)
  if (arg) {
    msg <- paste0(
      "Expected <", arg_name, "> to be FALSE or 0, but found ", arg
    )
    stop(msg)
  }
  TRUE
}

# checks that data frame contains a variable
# @param var_name the variable to be searched for
# @param data a data frame
# @param arg_name name of the data frame
check_in_data <- function(var_name, data, arg_name) {
  d_names <- colnames(data)
  ok <- (var_name %in% d_names)
  if (!ok) {
    str <- paste(d_names, collapse = ", ")
    msg <- paste0(
      "The variable '", var_name, "' not found in <", arg_name, ">!",
      " Found columns = {", str, "}."
    )
    stop(msg)
  }
  TRUE
}

# checks that argument is a \code{data.frame} and contains
# required variables (var_names)
check_df_with <- function(arg, var_names) {
  arg_name <- deparse(substitute(arg))
  c_data <- class(arg)
  if (c_data != "data.frame") {
    msg <- paste0("<", arg_name, "> must be a data.frame! found = ", c_data)
    stop(msg)
  }
  var_names <- unique(var_names)
  for (name in var_names) check_in_data(name, arg, arg_name)
  TRUE
}

# checks that argument has values less than equal to
# given maximums, elementwise
check_all_leq <- function(arg, maximums) {
  check_lengths(arg, maximums)
  L <- length(arg)
  for (j in seq_len(L)) {
    a <- arg[j]
    b <- maximums[j]
    if (a > b) {
      arg_name <- deparse(substitute(arg))
      m_name <- deparse(substitute(maximums))
      msg <- paste0(
        "value of <", arg_name, "> is larger than value of <",
        m_name, "> at index ", j, "!"
      )
      stop(msg)
    }
  }
  TRUE
}

# checks that object has no names
check_not_named <- function(arg) {
  arg_name <- deparse(substitute(arg))
  nams <- names(arg)
  if (!is.null(nams)) {
    str <- paste(nams, collapse = ", ")
    msg <- paste0("<", arg_name, "> should not have names! found = {", str, "}")
    stop(msg)
  }
  TRUE
}

# checks that object has names
check_named <- function(arg) {
  arg_name <- deparse(substitute(arg))
  nams <- names(arg)
  if (length(nams) < 1) {
    msg <- paste0("<", arg_name, "> must have names!")
    stop(msg)
  }
  TRUE
}

# checks that argument has expected number of dimensions
check_dim <- function(arg, D) {
  arg_name <- deparse(substitute(arg))
  L <- length(dim(arg))
  if (L != D) {
    msg <- paste0(
      "number of dimensions of <", arg_name, "> must be ", D,
      "! found = ", L
    )
    stop(msg)
  }
  TRUE
}

# checks that vector has only unique values
check_unique <- function(arg) {
  arg_name <- deparse(substitute(arg))
  L1 <- length(arg)
  L2 <- length(unique(arg))
  if (L1 != L2) stop("<", arg_name, "> must have only unique elements!")
  TRUE
}

# Check if argument is one of the allowed options
#
# @description Replacement for the \code{base} R function \code{match.arg}.
# Gives more informative errors and requires exact match. Should only be used
# as a subroutine of other functions and never directly.
#
# @param arg The given argument. Cannot be a list.
# @param allowed Allowed arguments. Must be have length at least 2.
# @return Return the index of \code{arg} in \code{allowed} or
# throws an error if argument is not valid.
# @family argument checks
check_allowed <- function(arg, allowed) {

  # Get names of given arguments and caller function
  arg_name <- deparse(substitute(arg))
  allowed_name <- deparse(substitute(allowed))
  caller_name <- deparse(sys.call(-1))
  msg <- paste0("Error in ", caller_name, ": ")

  # Check that 'arg' is not a list or vector
  if (length(arg) != 1) {
    msg <- paste0(msg, "Length of <", arg_name, "> must be one!")
    stop(msg)
  }

  # Check that there are at least two options
  if (length(allowed) < 2) {
    msg <- paste0(msg, "Length of <", allowed_name, "> must be at least 2!")
    stop(msg)
  }

  # Check that arg matches an allowed argument exactly
  idx <- which(allowed == arg)
  arg_msg <- paste0(
    "The given value '", arg,
    "' for argument <", arg_name, "> is invalid.\n"
  )

  L <- length(idx)
  if (L == 0) {
    msg <- paste0(msg, arg_msg, "Allowed values are {")
    opts <- paste0(allowed, collapse = ", ")
    msg <- paste0(msg, opts, "}.")
    stop(msg)
  } else if (L > 1) {
    msg <- paste0(msg, arg_msg, "It matches multiple allowed values.")
    stop(msg)
  }

  return(idx)
}
