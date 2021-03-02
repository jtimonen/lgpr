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
  check_not_null(var_name)
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

#' Return first list element if the list has length is one
#'
#' @param x a list
#' @return the original list or just its first element if its length is one
#' @family list utilities
simplify_list <- function(x) {
  check_type(x, "list")
  L <- length(x)
  if (L == 1) x <- x[[1]]
  return(x)
}
