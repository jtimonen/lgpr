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
  } else {
    check_dim(v, 0)
    n <- length(v)
    out <- array(v, dim = c(1, n))
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
  check_dim(x, 2)
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
  D_in <- dim(x)
  L <- length(D_in)
  if (L != 3) {
    msg <- paste0("Number of dimensions in <x> must be 3! found = ", L)
    stop(msg)
  }
  if (D_in[2] < 1) {
    msg <- paste0("2nd dimension of <x> must be at least 1! found = ", D_in[2])
    stop(msg)
  }
  D_out <- c(D_in[1] * D_in[2], D_in[3])
  out <- array(0, dim = D_out)
  for (j in seq_len(D_in[2])) {
    i1 <- (j - 1) * D_in[1] + 1
    i2 <- j * D_in[1]
    out[i1:i2, ] <- x[, j, ]
  }
  return(out)
}

#' Array transformation utility
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param L an integer
#' @return a list of length \code{L}, where each element is an array of
#' shape \code{h} x (\code{m}/\code{L}), where \code{h = nrow(x)}
#' @family array utilities
array_to_arraylist <- function(x, L) {
  m <- dim(x)[2]
  if (m %% L) stop("Second dimension of <x> must be a multiple of <L>!")
  out <- list()
  for (k in seq_len(L)) {
    inds <- seq(k, m, by = L)
    out[[k]] <- ensure_2dim(x[, inds])
  }
  return(out)
}

#' Add a vector to each column of an array
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param v a vector of length \code{n}
#' @return an array of shape \code{n} x \code{m}
#' @family array utilities
add_to_columns <- function(x, v) {
  check_dim(x, 2)
  check_length(v, nrow(x))
  D <- ncol(x)
  for (j in seq_len(D)) x[, j] <- x[, j] + v
  return(x)
}

#' Apply possible reduction to rows of an array
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param reduce a function or \code{NULL}
#' @return an array of shape \code{1} x \code{m} if \code{reduce} is not
#' \code{NULL}, otherwise the original array \code{x}
#' @family array utilities
apply_reduce <- function(x, reduce) {
  if (!is.null(reduce)) {
    check_type(reduce, "function")
    x <- apply(x, 2, reduce)
    x <- ensure_2dim(x)
  }
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
  out <- if (S < L) x else NULL
  return(out)
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

#' Helper function for creating plot data frames from prediction objects
#'
#' @param x an array with \code{num_draws} rows
#' @return a data frame with \code{num_draws} columns
#' @family array utilities
make_draw_df <- function(x) {
  num_draws <- nrow(x)
  x <- data.frame(t(x))
  colnames(x) <- paste("draw", seq_len(num_draws), sep = "_")
  rownames(x) <- NULL
  return(x)
}

#' Format a 3-dimensional array into a list of 2-dimensional arrays
#'
#' @param x an array
#' @return a list of arrays
#' @family array utilities
arr3_to_list <- function(x) {
  out <- list()
  d <- dim(x)[1]
  for (j in seq_len(d)) {
    out[[j]] <- arr3_select(x, j)
  }
  return(out)
}

#' Select an element along the first dimension of an array, and only drop
#' the first dimension
#'
#' @param x an array of shape\code{d1} x \code{d2} x \code{d3}
#' @param idx the index to select
#' @return an array of shape \code{d2} x \code{d3}
#' @family array utilities
arr3_select <- function(x, idx) {
  DIM <- dim(x)
  array(x[idx, , ], dim = DIM[2:3])
}
