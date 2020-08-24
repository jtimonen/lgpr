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
#' @param x a listof length \code{m} where each element is a vector of
#' length code \code{n}
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
  index <- argument_check(likelihood, allowed)
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

#' Get covariate names for a model
#'
#' @param model an object of class \linkS4class{lgpmodel}
#' @param type what types of covariates to include. Must be either
#' "all", "categorical" or "continuous".
#' @return A string where names are separated by a comma.
covariate_names <- function(model, type = "all") {
  allowed <- c("continuous", "categorical", "all")
  nam1 <- model@var_names$x_cont
  nam2 <- model@var_names$x_cat
  nam3 <- c(nam1, nam2)
  names <- list(nam1, nam2, nam3)
  idx <- argument_check(arg = type, allowed = allowed)
  out <- names[[idx]]
  paste(out, collapse = ", ")
}
