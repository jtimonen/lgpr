#' Check if argument has correct type.
#'
#' @param object Any object.
#' @param allowed Allowed class names.
#' @return Returns \code{TRUE} if the object has an allowed type.
check_type <- function(object, allowed) {
  type <- class(object)
  ok <- (type %in% allowed)
  if (!ok) {
    arg_name <- deparse(substitute(object))
    str <- paste(allowed, collapse = ", ")
    msg <- paste0(
      arg_name, " has invalid type '", type,
      "'. Allowed types are {", str, "}."
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

#' Check if a value is numeric and throw error if not
#'
#' @param arg A single number
#' @param arg_name Argument name
#' @param require_positive Should this also check that \code{arg} is positive?
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

#' Check if argument has correct length.
#'
#' @param object Any object.
#' @param expected Expected length.
#' @return Returns \code{TRUE} if the object has an allowed length.
check_length <- function(object, expected) {
  L <- length(object)
  ok <- (L == expected)
  if (!ok) {
    arg_name <- deparse(substitute(object))
    msg <- paste0(
      arg_name, " has length ", L, ", but its length should be ",
      expected, "!"
    )
    stop(msg)
  }
  return(TRUE)
}
