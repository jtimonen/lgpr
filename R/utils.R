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
      "Element with name '", var_name,
      "' not found in <", obj_name, ">!"
    )
    msg <- paste0(msg, " Found elements: {", valid, "}")
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

# -------------------------------------------------------------------------

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

# -------------------------------------------------------------------------

#' Helper function for generic functions
#'
#' @description Helper function for generic functions that work on
#' both of \linkS4class{lgpmodel} and \linkS4class{lgpfit} class objects.
#' @param object an object of class \linkS4class{lgpmodel} or
#' \linkS4class{lgpfit}
#' @return an object of class \linkS4class{lgpmodel}
object_to_model <- function(object) {
  allowed <- c("lgpmodel", "lgpfit")
  check_type(object, allowed)
  if (class(object) == "lgpfit") {
    out <- object@model
  } else {
    out <- object
  }
  return(out)
}

#' Infer observed effect times from a data frame
#'
#' @param data a data frame
#' @param age_variable age variable name
#' @param disage_variable disease-related age variable name
#' @param id_variable id variable name
#' @return a data frame
get_teff_obs <- function(data,
                         age_variable,
                         disage_variable,
                         id_variable) {
  check_type(data, "data.frame")
  fac <- dollar(data, id_variable)
  t1 <- dollar(data, age_variable)
  t2 <- dollar(data, disage_variable)
  uid <- unique(fac)
  L <- length(uid)
  teff <- rep(0, L)
  nam <- rep(0, L)
  for (j in seq_len(L)) {
    bbb <- as.numeric(fac == uid[j]) + as.numeric(t2 == 0)
    i0 <- which(bbb == 2)
    L <- length(i0)
    if (L > 0) {
      idx <- i0[1]
      teff[j] <- t1[idx]
    } else {
      teff[j] <- NaN
    }
    nam[j] <- uid[j]
  }
  df <- data.frame(nam, teff)
  colnames(df) <- c(id_variable, age_variable)
  return(df)
}

#' Integer encoding of likelihod function names
#'
#' @description
#' \itemize{
#'   \item \code{likelihood_as_int} converts likelihood name to Stan encoding
#'   \item \code{likelihood_as_str} converts the Stan likelihood encoding
#'   to a string
#'   \item \code{likelihood_list} returns the available likelihood names
#' }
#' @param likelihood a string
#' @param index an integer
#' @name likelihood_encoding

#' @rdname likelihood_encoding
likelihood_list <- function() {
  c("gaussian", "poisson", "nb", "binomial", "bb")
}

#' @rdname likelihood_encoding
likelihood_as_str <- function(index) {
  names <- likelihood_list()
  name <- names[index]
  if (is.na(name)) {
    stop("invalid likelihood index")
  }
  name
}

#' @rdname likelihood_encoding
likelihood_as_int <- function(likelihood) {
  likelihood <- tolower(likelihood)
  allowed <- likelihood_list()
  index <- check_allowed(likelihood, allowed)
  return(index)
}

#' Link functions and their inverses
#'
#' @param x the input
#' @param likelihood name of the likelihood model
#' @returns transformed input
#' @name link
NULL

#' @rdname link
link <- function(x, likelihood) {
  allowed <- likelihood_list()
  check_allowed(likelihood, allowed)
  if (likelihood %in% c("poisson", "nb")) {
    x <- log(x)
  } else if (likelihood %in% c("binomial", "bb")) {
    x <- log(x) - log(1 - x)
  }
  return(x)
}

#' @rdname link
link_inv <- function(x, likelihood) {
  allowed <- likelihood_list()
  check_allowed(likelihood, allowed)
  if (likelihood %in% c("poisson", "nb")) {
    x <- exp(x)
  } else if (likelihood %in% c("binomial", "bb")) {
    x <- 1 / (1 + exp(-x))
  }
  return(x)
}

#' Get lgpr version description
#'
#' @export
#' @return package description
get_pkg_description <- function() {
  lgprLib <- dirname(system.file(package = "lgpr"))
  descr <- suppressWarnings(utils::packageDescription("lgpr",
    lib.loc = lgprLib
  ))
  return(descr)
}

#' Get a stan model of the package
#'
#' @export
#' @param name name of the model
#' @return an object of class stanmodel
get_stan_model <- function(name = "lgp") {
  dollar(stanmodels, name)
}

#' Ensure vector has expected length
#'
#' @param len the expected length
#' @param v original vector or just one value that is replicated
#' @return a vector of length \code{len}
ensure_len <- function(v, len) {
  v_name <- deparse(substitute(v))
  L <- length(v)
  if (L == 1) {
    v <- rep(v, len)
  } else if (L != len) {
    msg <- paste0(
      "length of <", v_name, "> was expected to be 1 or ", len,
      ", but found length ", L
    )
    stop(msg)
  }
  return(v)
}
