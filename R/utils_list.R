#' A safer alternative for the dollar operator
#'
#' @description Requires exact match and throws error if \code{var_name} is
#' not in \code{names(object)}.
#' @param object a list or data frame
#' @param var_name name of the variable to access
#' @return Returns \code{object[[var_name, exact = TRUE]]} if variable
#' exists.
#' @family list utilities
dollar <- function(object, var_name) {
  obj_name <- deparse(substitute(object))
  nams <- names(object)
  if (!(var_name %in% nams)) {
    valid <- paste(nams, collapse = ", ")
    msg <- paste0(
      "Variable with name '", var_name,
      "' not found in <", obj_name, ">!"
    )
    msg <- paste0(msg, " Found: {", valid, "}")
    stop(msg)
  }
  object[[var_name, exact = TRUE]]
}

#' Printing a list in a more compact format
#'
#' @param input a named list
#' @family list utilities
print_list <- function(input) {

  # Helper function
  is_printable <- function(object) {
    d <- dim(object)
    if (is.null(d)) {
      return(TRUE)
    } else if (prod(d) == 0) {
      return(FALSE)
    }
    TRUE
  }

  # Choose which fields to print
  nam <- names(input)
  printed <- c()
  skipped <- c()
  for (name in nam) {
    f <- input[[name]]
    if (is_printable(f)) {
      printed <- c(printed, name)
    } else {
      skipped <- c(skipped, name)
    }
  }

  # Print
  print(input[printed])
  str <- paste(skipped, collapse = ", ")
  msg <- paste0(
    "Did not print fields with at least one zero dimension:\n    ",
    str, "\n"
  )
  cat(msg)
  invisible(input)
}

#' Zip two lists into a list of lists of length 2
#'
#' @param a a list of length \code{L}
#' @param b a list of length \code{L}
#' @return a list of length \code{L}
#' @family list utilities
zip_lists <- function(a, b) {
  check_type(a, "list")
  check_type(b, "list")
  check_lengths(a, b)
  a_name <- deparse(substitute(a))
  b_name <- deparse(substitute(b))
  out <- list()
  L <- length(a)
  for (j in seq_len(L)) {
    list_j <- list(a[[j]], b[[j]])
    names(list_j) <- c(a_name, b_name)
    out[[j]] <- list_j
  }
  return(out)
}

#' Wrap list into a list of length 1 if the original list is named
#'
#' @param x a list
#' @return a list with no names
#' @family list utilities
list_if_named <- function(x) {
  check_type(x, "list")
  is_named <- !is.null(names(x))
  if (is_named) x <- list(x)
  return(x)
}

#' List elements to matrix rows
#'
#' @param x a list of length \code{m} where each element is a vector of
#' length \code{n}
#' @param n length of each vector
#' @return a matrix with shape \code{m} x \code{n}
#' @family list utilities
list_to_matrix <- function(x, n) {
  m <- length(x)
  A <- array(0, dim = c(m, n))
  for (i in seq_len(m)) {
    A[i, ] <- x[[i]]
  }
  as.matrix(A)
}

#' Matrix rows to a list
#'
#' @param x a matrix or array with \code{m} rows and \code{n} columns
#' @return an unnamed list of length \code{m} where each element is a
#' vector of length \code{n}
#' @family list utilities
matrix_to_list <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- x[i, ]
  }
  return(L)
}
