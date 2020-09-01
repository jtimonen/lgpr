#' Argument validation
#'
#' @name checks
#' @param arg An argument to check.
#' @description
#' \itemize{
#'   \item \code{check_type} checks if argument has correct class
#'   \item \code{check_numeric} checks if argument is numeric
#'   \item \code{check_not_null} checks if argument is not null
#'   \item \code{check_false} checks if argument is false or zero
#'   \item \code{check_length} checks if argument has correct length
#'   \item \code{check_lengths} checks if two argument have equal length
#' }
#' @return \code{TRUE} if the check passes.
#' @family argument checks
NULL

#' @rdname checks
#' @param allowed Allowed class names.
check_type <- function(arg, allowed) {
  check_length(class(arg), 1)
  type <- class(arg)
  ok <- (type %in% allowed)
  if (!ok) {
    arg_name <- deparse(substitute(arg))
    str <- paste(allowed, collapse = ", ")
    msg <- paste0(
      arg_name, " has invalid type '", type,
      "'. Allowed types are {", str, "}."
    )
    stop(msg)
  }
  return(TRUE)
}

#' @rdname checks
#' @param arg_name Argument name.
#' @param require_positive Should it also be checked that \code{arg} is
#' positive?
check_numeric <- function(arg, arg_name = NULL, require_positive = FALSE) {
  if (is.null(arg_name)) {
    arg_name <- deparse(substitute(arg))
  }
  if (!is.numeric(arg)) {
    msg <- paste0("<", arg_name, "> must be numeric! found = ", arg)
    stop(msg)
  } else if (arg <= 0) {
    if (require_positive) {
      msg <- paste0("<", arg_name, "> must be positive! found = ", arg)
      stop(msg)
    }
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
check_not_null <- function(arg) {
  arg_name <- deparse(substitute(arg))
  if (is.null(arg)) {
    stop(arg_name, " should not be NULL!")
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
