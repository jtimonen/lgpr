#' Argument validation
#'
#' @name checks
#' @param arg An argument to check.
#' @description
#' \itemize{
#'   \item \code{check_type} checks if argument has correct class
#'   \item \code{check_function} checks if argument is a function
#'   \item \code{check_numeric} checks if argument is numeric
#'   \item \code{check_positive} checks if argument is positive
#'   \item \code{check_non_negative} checks if argument is non-negative
#'   \item \code{check_positive_all} checks if argument has only positive
#'   values
#'   \item \code{check_non_negative_all} checks if argument has only
#'   non-negative values
#'   \item \code{check_integer_all} checks if argument has only
#'   integer values
#'   \item \code{check_interval} checks if argument is inside a given interval
#'   \item \code{check_not_null} checks if argument is not null
#'   \item \code{check_false} checks if argument is false or zero
#'   \item \code{check_length} checks if argument has given length
#'   \item \code{check_length_1_or} checks if argument has given length or
#'   length 1
#'   \item \code{check_length_geq} checks if argument has at least a
#'   given length
#'   \item \code{check_lengths} checks if two arguments have equal length
#'   \item \code{check_in_data} checks that variable exists in a data frame
#'   \item \code{check_all_leq} checks that argument has values less than equal
#'   to given maximums, elementwise
#'   \item \code{check_not_named} checks that \code{names(arg)} is NULL
#'   \item \code{check_named} checks that length of \code{names(arg)} is
#'   greater than zero
#'   \item \code{check_dim} checks that \code{arg} has expected number of
#'   dimensions
#'   \item \code{check_null} checks that \code{arg} is NULL
#' }
#' @return \code{TRUE} if the check passes.
#' @family argument checks
NULL

#' @rdname checks
#' @param allowed Allowed class names.
check_type <- function(arg, allowed) {
  arg_name <- deparse(substitute(arg))

  # If allowed is 'function'
  if (length(allowed) == 1) {
    if (allowed == "function") {
      check_function(arg, arg_name)
      return(TRUE)
    }
  }

  # Check type otherwise
  type <- class(arg)
  ok <- any(type %in% allowed)
  if (!ok) {
    str1 <- paste(allowed, collapse = ", ")
    str2 <- paste(type, collapse = ", ")
    msg <- paste0(
      "Wrong argument type: class(", arg_name, ") must contain one of {",
      str1, "}. Found = {", str2, "}."
    )
    stop(msg)
  }

  return(TRUE)
}

#' @rdname checks
#' @param arg_name argument name
check_function <- function(arg, arg_name) {
  ok <- is.function(arg)
  if (!ok) {
    msg <- paste0(
      "'", arg_name, "' must be a function, but ",
      "is.function(", arg_name, ") returned FALSE."
    )
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
check_numeric <- function(arg) {
  if (!is.numeric(arg)) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must be numeric! found = ", arg)
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
check_positive <- function(arg) {
  check_numeric(arg)
  if (arg <= 0) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must be positive! found = ", arg)
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
check_non_negative <- function(arg) {
  check_numeric(arg)
  if (arg < 0) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must be non-negative! found = ", arg)
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
check_positive_all <- function(arg) {
  check_numeric(arg)
  if (any(arg <= 0)) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must have only positive values")
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
check_non_negative_all <- function(arg) {
  check_numeric(arg)
  if (any(arg < 0)) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must have only non-negative values")
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
check_integer_all <- function(arg) {
  check_numeric(arg)
  if (sum(round(arg) != arg) > 0) {
    arg_name <- deparse(substitute(arg))
    msg <- paste0("<", arg_name, "> must have only integer values")
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
#' @param lower interval minimum
#' @param upper interval maximum
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
  return(TRUE)
}

#' @rdname checks
check_not_null <- function(arg) {
  arg_name <- deparse(substitute(arg))
  if (is.null(arg)) {
    stop(arg_name, " should not be NULL!")
  }
  return(TRUE)
}

#' @rdname checks
#' @param reason explanation why the argument should be NULL
check_null <- function(arg, reason = NULL) {
  arg_name <- deparse(substitute(arg))
  if (!is.null(arg)) {
    msg <- paste0("<", arg_name, "> should be NULL!")
    if (!is.null(reason)) {
      msg <- paste0(msg, " Reason: ", reason)
    }
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
#' @param expected Expected length.
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
  return(TRUE)
}

#' @rdname checks
#' @param expected Expected length.
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
  return(TRUE)
}

#' @rdname checks
#' @param min_length Minimum allowed length.
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
  return(TRUE)
}

#' @param a An argument to check.
#' @param b An argument to check.
#' @rdname checks
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
  return(TRUE)
}

#' @rdname checks
check_false <- function(arg) {
  arg_name <- deparse(substitute(arg))
  check_not_null(arg)
  if (arg) {
    msg <- paste0(
      "Expected <", arg_name, "> to be FALSE or 0, but found ", arg
    )
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
#' @param var_name the variable to be searched for
#' @param data a data frame
check_in_data <- function(var_name, data) {
  d_names <- colnames(data)
  ok <- (var_name %in% d_names)
  if (!ok) {
    str <- paste(d_names, collapse = ", ")
    msg <- paste0(
      "The variable '", var_name, "' not found in <data>!",
      " Found data columns = {", str, "}."
    )
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
#' @param maximums maximum allowed values for \code{arg}
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

#' @rdname checks
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

#' @rdname checks
check_named <- function(arg) {
  arg_name <- deparse(substitute(arg))
  nams <- names(arg)
  if (length(nams) < 1) {
    msg <- paste0("<", arg_name, "> must have names!")
    stop(msg)
  }
  TRUE
}

#' @rdname checks
#' @param D expected number of dimensions
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

#' Check if argument is one of the allowed options
#'
#' @description Replacement for the \code{base} R function \code{match.arg}.
#' Gives more informative errors and requires exact match. Should only be used
#' as a subroutine of other functions and never directly.
#'
#' @param arg The given argument. Cannot be a list.
#' @param allowed Allowed arguments. Must be have length at least 2.
#' @return Return the index of \code{arg} in \code{allowed} or
#' throws an error if argument is not valid.
#' @family argument checks
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
