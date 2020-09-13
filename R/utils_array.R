#' Repeat a vector as a rows of an array
#'
#' @description Throws an error if \code{v} is \code{NULL}.
#' @param v a vector of length \code{m}
#' @param n number of times to repeat
#' @return returns an array of size \code{n} x \code{m}
#' @family array utilities
repvec <- function(v, n) {
  check_not_null(v)
  m <- length(v)
  A <- matrix(rep(v, n), n, m, byrow = TRUE)
  return(as.array(A))
}

#' Ensure that v is 2-dimensional
#'
#' @param v a vector of length \code{m} or an array of size \code{n} x \code{m}
#' @return returns a 2-dimensional array
#' @family array utilities
ensure_2dim <- function(v) {
  check_not_null(v)
  L <- length(dim(v))
  if (L == 2) {
    out <- v
  } else if (L < 2) {
    n <- length(v)
    out <- array(v, dim = c(1, n))
  } else {
    stop("<v> must have one or two dimensions!")
  }
  return(out)
}

#' Compute row variances for a 2-dimensional array
#'
#' @param x an array of size \code{n} x \code{m}
#' @return returns a vector
#' @family array utilities
row_vars <- function(x) {
  check_not_null(x)
  L <- length(dim(x))
  if (L != 2) stop("<x> must be 2-dimensional")
  apply(x, 1, stats::var)
}

#' Normalize matrix or data frame rows so that they sum to 1
#'
#' @param x an array of size \code{n} x \code{m}
#' @return an array of size \code{n} x \code{m}
#' @family array utilities
normalize_rows <- function(x) {
  s <- rowSums(x)
  x / s
}

#' Select named row of an array
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param name name of the row
#' @return a vector of length \code{m}
#' @family array utilities
select_row <- function(x, name) {
  df <- data.frame(t(x))
  dollar(df, name)
}

#' Squeeze the second dimension of an array
#'
#' @param x an array of shape \code{n1} x \code{n2} x \code{n3}
#' @return an array of shape \code{n1*n2} x \code{n3}
#' @family array utilities
squeeze_second_dim <- function(x) {
  D <- dim(x)
  L <- length(D)
  if (L != 3) {
    msg <- paste0("Number of dimensions in <x> must be 3! found = ", L)
    stop(msg)
  }
  if (D[2] < 1) {
    msg <- paste0("Second dimension of <x> must be at least 1! found = ", D[2])
    stop(msg)
  }
  a <- c()
  for (j in seq_len(D[2])) {
    b <- x[, j, ]
    a <- rbind(a, b)
  }
  return(a)
}

#' Array reshaping
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param L an integer
#' @param draws see the \code{draws} argument of \code{\link{get_f}}
#' @return a list of length \code{L}, where each element is an array of
#' shape \code{h} x (\code{m}/\code{L}), where \code{h = length(draws)}
#' @family array utilities
array_to_arraylist <- function(x, L, draws) {
  m <- dim(x)[2]
  if (m %% L) {
    stop("Second dimension of <x> must be a multiple of <L>!")
  }
  out <- list()
  for (k in seq_len(L)) {
    inds <- seq(k, m, by = L)
    out[[k]] <- ensure_2dim(x[draws, inds])
  }
  return(out)
}

#' Add the sum of all arrays in a list to the list
#'
#' @param x  a list of arrays of shape \code{n} x \code{m}, with
#' length \code{L}
#' @return a list of arrays of shape \code{n} x \code{m},
#' with length \code{L + 1}
#' @family array utilities
add_sum_arraylist <- function(x) {
  L <- length(x)
  if (L == 0) stop("list has length 0!")
  s <- x[[1]]
  for (j in 2:L) {
    s <- s + x[[j]]
  }
  x[[L + 1]] <- s
  return(x)
}

#' Return NULL if vector contains only NaN values
#'
#' @param x a vector
#' @return \code{x} unchanged or \code{NULL}
#' @family array utilities
null_if_all_nan <- function(x) {
  L <- length(x)
  S <- sum(is.nan(x))
  if (S < L) {
    return(x)
  } else {
    return(NULL)
  }
}

#' Check that each row of array is identical
#'
#' @param rows rows
#' @return the indentical row
#' @family array utilities
reduce_rows <- function(rows) {
  nam <- rownames(rows)
  R <- dim(rows)[1]
  r1 <- as.numeric(rows[1, ])
  nam1 <- nam[1]
  n <- length(r1)
  if (R == 1) {
    return(r1)
  }
  for (j in 2:R) {
    rj <- as.numeric(rows[j, ])
    namj <- nam[j]
    s <- sum(rj == r1)
    if (s != n) {
      msg <- paste0(
        "For each term with uncrt() or heter() expressions, ",
        "NaNs of the continuous covariate must be on the same",
        " rows. Found discrepancy between ",
        nam1, " and ", namj, "."
      )
      stop(msg)
    }
  }
  return(r1)
}
