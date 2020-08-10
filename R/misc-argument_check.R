#' Check if argument is valid
#'
#' @description Replacement for the \code{base} R function \code{match.arg}.
#' Gives more informative errors and requires exact match. Should only be used
#' as a subroutine of other functions and never directly.
#'
#' @param arg The given argument. Cannot be a list.
#' @param allowed Allowed arguments. Must be have length at least 2.
#' @return Return the index of \code{arg} in \code{allowed} or
#' throws an error if argument is not valid.
argument_check <- function(arg, allowed) {

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
