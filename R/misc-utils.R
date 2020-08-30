#' Repeat a vector as a rows of an array
#'
#' @description Throws an error if \code{v} is \code{NULL}.
#' @param v a vector of length \code{m}
#' @param n number of times to repeat
#' @return returns an array of size \code{n} x \code{m}
repvec <- function(v, n) {
  if (is.null(v)) {
    stop("Error in repvec, input <v> is NULL.")
  }
  m <- length(v)
  A <- matrix(rep(v, n), n, m, byrow = TRUE)
  return(as.array(A))
}

#' Access list field and throw error if it is NULL
#'
#' @param fields list
#' @param name field name
#' @return Returns \code{fields[[name]]} unless it is \code{NULL}
field_must_exist <- function(fields, name) {
  desc <- fields[[name]]
  arg_name <- deparse(substitute(fields))
  if (is.null(desc)) {
    str <- paste(names(fields), collapse = ", ")
    msg <- paste0(
      "The list <", arg_name, "> must contain a field named ",
      name, "! Found = {", str, "}"
    )
    stop(msg)
  }
  return(desc)
}

#' Wrapper for rstan::get_stream
#'
#' @description See \code{\link[rstan]{get_stream}}.
#' @return an external pointer
get_stream <- function() {
  rstan::get_stream()
}

#' List elements to matrix rows
#'
#' @param x a list of length \code{m} where each element is a vector of
#' length \code{n}
#' @param n length of each vector
#' @return a matrix with shape \code{m} x \code{n}
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
matrix_to_list <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- x[i, ]
  }
  return(L)
}

#' Names that the list given as data to Stan should contain
#'
#' @return a character vector
stan_list_names <- function() {
  c(
    "is_verbose",
    "is_generated_skipped",
    "is_f_sampled",
    "is_likelihood_skipped",

    "num_obs",
    "num_cov_cont",
    "num_cov_cat",
    "num_comps",
    "num_ell",
    "num_ns",
    "num_heter",
    "num_uncrt",
    "num_cases",

    "obs_model",
    "components",
    "y_cont",
    "y_disc",
    "y_num_trials",
    "x_cat",
    "x_cat_num_levels",
    "x_cont",
    "x_cont_unnorm",
    "x_cont_mask",

    "c_hat",
    "delta",
    "vm_params",
    "idx_expand",

    "prior_alpha",
    "prior_ell",
    "prior_wrp",
    "prior_sigma",
    "prior_phi",
    "prior_teff",

    "hyper_alpha",
    "hyper_ell",
    "hyper_wrp",
    "hyper_phi",
    "hyper_sigma",
    "hyper_beta",
    "hyper_teff",

    "teff_obs",
    "teff_lb",
    "teff_ub"
  )
}

#' List of allowed observation models
#'
#' @return a list of likelihood names
likelihood_list <- function() {
  c("gaussian", "poisson", "nb", "binomial")
}

#' Convert the Stan likelihood encoding to a string
#'
#' @param index an integer
#' @return a string
likelihood_as_str <- function(index) {
  names <- likelihood_list()
  name <- names[index]
  if (is.na(name)) {
    stop("invalid likelihood index")
  }
  name
}

#' Convert likelihood string to Stan encoding
#'
#' @param likelihood a string
#' @return an integer
likelihood_as_int <- function(likelihood) {
  likelihood <- tolower(likelihood)
  allowed <- likelihood_list()
  index <- check_allowed(likelihood, allowed)
  return(index)
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

#' Get main stan model of the package
#'
#' @export
#' @return an object of class stanmodel
get_stan_model <- function() {
  return(stanmodels[["lgp"]])
}

#' Default variance masking function parameters are defined here
#'
#' @return two numbers
default_vm_params <- function() {
  c(0.025, 1)
}


#' Return NULL if vector contains only NaN values
#'
#' @param x a vector
#' @return \code{x} unchanged or \code{NULL}
null_if_all_nan <- function(x) {
  L <- length(x)
  S <- sum(is.nan(x))
  if (S < L) {
    return(x)
  } else {
    return(NULL)
  }
}

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

#' Class info for show methods
#'
#' @param class_name class name
#' @return a string
class_info <- function(class_name) {
  str <- paste0(
    "An object of class ", class_name, ". See ?",
    class_name, " for more info and related methods/functions.\n"
  )
  return(str)
}

#' Squeeze the second dimension of an array
#'
#' @param x an array of shape \code{n1} x \code{n2} x \code{n3}
#' @return an array of shape \code{n1*n2} x \code{n3}
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
#' @param draws see the \code{draws} argument of \code{\link{get_posterior_f}}
#' @return a list of length \code{L}, where each element is an array of
#' shape \code{h} x (\code{m}/\code{L}), where \code{h = length(draws)}
array_to_arraylist <- function(x, L, draws) {
  m <- dim(x)[2]
  if (m %% L) {
    stop("Second dimension of <x> must be a multiple of <L>!")
  }
  out <- list()
  for (k in seq_len(L)) {
    inds <- seq(k, m, by = L)
    out[[k]] <- x[draws, inds]
  }
  return(out)
}

#' Is an object printable
#'
#' @param object an object
is_printable <- function(object) {
  d <- dim(object)
  if (is.null(d)) {
    return(TRUE)
  } else if (prod(d) == 0) {
    return(FALSE)
  }
  TRUE
}

#' Print a list in a more compact format
#'
#' @param input a named list
print_list <- function(input) {
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
zip_lists <- function(a, b) {
  check_type(a, "list")
  check_type(b, "list")
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
  out <- list()
  for (j in seq_len(L1)) {
    list_j <- list(a[[j]], b[[j]])
    names(list_j) <- c(a_name, b_name)
    out[[j]] <- list_j
  }
  return(out)
}
