#' Parse the covariates and model components from given data and formula
#'
#' @inheritParams parse_response
#' @return parsed input to stan and covariate scaling
parse_covs_and_comps <- function(data, model_formula) {

  # Check that all covariates exist in data
  x_names <- rhs_variables(model_formula@terms)
  x_names <- unique(x_names)
  for (name in x_names) {
    check_in_data(name, data)
  }

  # Create the inputs to Stan
  covariates <- stan_data_covariates(data, x_names)
  components <- stan_data_components(model_formula, covariates)
  to_stan <- c(covariates$to_stan, components$to_stan)

  # Return
  list(
    to_stan = to_stan,
    x_cont_scalings = covariates$x_cont_scalings,
    x_cat_levels = covariates$x_cat_levels
  )
}

#' Create covariate data for Stan input
#'
#' @description Creates the following Stan data input list fields:
#' \itemize{
#'   \item \code{num_cov_cat}
#'   \item \code{num_cov_cont}
#'   \item \code{x_cat}
#'   \item \code{x_cat_num_levels}
#'   \item \code{x_cont}
#'   \item \code{x_cont_unnorm}
#'   \item \code{x_cont_mask}
#' }
#' @param data a data frame
#' @param x_names unique covariate names
#' @return a named list with fields
#' \itemize{
#'   \item \code{to_stan}: a list of stan data
#'   \item \code{x_cont_scaling}: normalization function and inverse for each
#'   continuous covariate
#'   \item \code{x_cat_levels}: names of the levels of each categorical
#'   covariate before conversion from factor to numeric
#' }
stan_data_covariates <- function(data, x_names) {
  num_obs <- dim(data)[1]

  x_cont <- list()
  x_cont_mask <- list()
  x_cont_unnorm <- list()
  x_cont_scalings <- list()
  x_cont_names <- c()

  x_cat <- list()
  x_cat_levels <- list()
  x_cat_num_levels <- c()
  x_cat_names <- c()

  num_cat <- 0
  num_cont <- 0

  for (name in x_names) {
    X_RAW <- data[[name]]
    c_x <- class(X_RAW)
    if (c_x == "factor") {

      # A categorical covariate
      num_cat <- num_cat + 1
      n_na <- sum(is.na(X_RAW))
      if (n_na > 0) {
        msg <- paste0(n_na, " missing values for factor '", name, "'!")
        stop(msg)
      }
      x_cat[[num_cat]] <- as.numeric(X_RAW)
      x_cat_num_levels[num_cat] <- length(levels(X_RAW))
      x_cat_levels[[num_cat]] <- levels(X_RAW)
      x_cat_names[num_cat] <- name
    } else if (c_x == "numeric") {

      # A continuous covariate
      num_cont <- num_cont + 1
      is_na <- is.na(X_RAW)
      x_cont_mask[[num_cont]] <- as.numeric(is_na)
      X_NONAN <- X_RAW
      X_NONAN[is_na] <- 0
      normalizer <- create_scaling(X_NONAN, name)
      x_cont_scalings[[num_cont]] <- normalizer
      x_cont[[num_cont]] <- normalizer@fun(X_NONAN)
      x_cont_unnorm[[num_cont]] <- X_NONAN
      x_cont_names[num_cont] <- name
    } else {
      msg <- paste0(
        "Covariate '", name, "' has invalid type '", c_x,
        "'! Must be one of {'factor', 'numeric'}"
      )
      stop(msg)
    }
  }

  # Convert lists to matrices
  x_cat <- list_to_matrix(x_cat, num_obs)
  x_cont <- list_to_matrix(x_cont, num_obs)
  x_cont_unnorm <- list_to_matrix(x_cont_unnorm, num_obs)
  x_cont_mask <- list_to_matrix(x_cont_mask, num_obs)

  # Name lists and matrix rows
  names(x_cont_scalings) <- x_cont_names
  names(x_cat_levels) <- x_cat_names
  rownames(x_cat) <- x_cat_names
  rownames(x_cont) <- x_cont_names
  rownames(x_cont_unnorm) <- x_cont_names
  rownames(x_cont_mask) <- x_cont_names

  if (!is.null(x_cat_num_levels)) {
    x_cat_num_levels <- array(x_cat_num_levels, dim = c(num_cat))
  } else {
    x_cat_num_levels <- array(0, dim = c(0))
  }


  # Create Stan data
  to_stan <- list(
    num_cov_cont = num_cont,
    num_cov_cat = num_cat,
    x_cat = x_cat,
    x_cat_num_levels = x_cat_num_levels,
    x_cont = x_cont,
    x_cont_unnorm = x_cont_unnorm,
    x_cont_mask = x_cont_mask
  )

  # Return
  list(
    to_stan = to_stan,
    x_cont_scalings = x_cont_scalings,
    x_cat_levels = x_cat_levels
  )
}

#' Create model components data for Stan input
#'
#' @param model_formula an object of class \linkS4class{lgpformula}
#' @param covariates a list returned by \code{\link{stan_data_covariates}}
#' @return a named list with fields
#' \itemize{
#'   \item \code{to_stan}: a list of stan data
#' }
stan_data_components <- function(model_formula, covariates) {
  terms <- model_formula@terms@summands
  J <- length(terms)

  # Create the components integer array
  comps <- array(0, dim = c(J, 9))
  for (j in seq_len(J)) {
    comps[j, ] <- term_to_numeric(terms[[j]], covariates)
  }
  colnames(comps) <- c(
    "type", "ker", "unused",
    "het", "ns", "vm",
    "unc", "icat", "icont"
  )
  rownames(comps) <- term_names(model_formula@terms)
  components <- as.matrix(comps)

  # Create idx_expand, num_cases and vm_params
  x_cat <- covariates$to_stan$x_cat
  x_cont_mask <- covariates$to_stan$x_cont_mask
  idx_expand <- create_idx_expand(components, x_cat, x_cont_mask)
  num_cases <- length(unique(idx_expand[idx_expand != 0]))
  num_ns <- sum(components[, 5] != 0)
  num_vm <- sum(components[, 6] != 0)
  VM <- default_vm_params()
  vm_params <- matrix(rep(VM, num_ns), num_ns, 2, byrow = TRUE)

  to_stan <- list(
    components = components,
    idx_expand = idx_expand,
    num_cases = num_cases,
    num_ell = sum(components[, 1] != 0),
    num_heter = sum(components[, 4] != 0),
    num_ns = num_ns,
    num_vm = num_vm,
    num_uncrt = sum(components[, 7] != 0),
    num_comps = dim(components)[1],
    vm_params = vm_params
  )

  # Return
  list(
    to_stan = to_stan
  )
}


#' Map a list of terms to their "names"
#'
#' @param rhs an object of class \linkS4class{lgprhs}
#' @return a character vector
term_names <- function(rhs) {
  terms <- rhs@summands
  J <- length(terms)
  names <- c()
  for (j in seq_len(J)) {
    term <- terms[[j]]
    names <- c(names, as.character(term))
  }
  return(names)
}


#' An lgpterm to numeric representation for Stan
#'
#' @param term an object of class \linkS4class{lgpterm}
#' @param covariates a list returned by \code{\link{stan_data_covariates}}
#' @return a vector of 9 integers
term_to_numeric <- function(term, covariates) {
  out <- rep(0, 9)

  # Check formula validity
  parsed <- check_term_factors(term)

  # Check component type
  is_gp <- !is.null(parsed$gp_kernel)
  is_cat <- !is.null(parsed$cat_kernel)
  if (!is_gp) {
    ctype <- 0
  } else {
    ctype <- if (is_cat) 2 else 1
  }
  out[1] <- ctype

  # Check kernel type
  if (is_cat) {
    kernels <- c("zerosum", "categ")
    idx <- check_allowed(parsed$cat_kernel, allowed = kernels)
    ktype <- idx - 1
  } else {
    ktype <- 0
  }
  out[2] <- ktype

  # Check nonstationary options
  gpk <- parsed$gp_kernel
  if (!is.null(gpk)) {
    is_warped <- parsed$gp_kernel %in% c("gp_warp", "gp_warp_vm")
    is_vm <- parsed$gp_kernel == "gp_warp_vm"
  } else {
    is_warped <- FALSE
    is_vm <- FALSE
  }

  out[5] <- as.numeric(is_warped)
  out[6] <- as.numeric(is_vm)

  # Check covariate types
  cidx <- check_term_covariates(covariates, parsed)
  out[4] <- cidx[1]
  out[7] <- cidx[2]
  out[8] <- cidx[3]
  out[9] <- cidx[4]

  return(out)
}

#' Helper for converting an lgpterm to numeric representation for Stan
#'
#' @param covariates a list returned by \code{\link{stan_data_covariates}}
#' @param pf a list returned by \code{\link{check_term_factors}}
#' @return two integers
check_term_covariates <- function(covariates, pf) {
  cat_names <- rownames(covariates$to_stan$x_cat)
  cont_names <- rownames(covariates$to_stan$x_cont)
  nams <- list(categorical_names = cat_names, continuous_names = cont_names)

  check_type <- function(cov_name, allowed, type, fun) {
    idx <- which(allowed == cov_name)
    if (length(idx) < 1) {
      msg <- paste0(
        "The argument for <", fun, "> must be a name of a ",
        type, " covariate. Found = ", cov_name, "."
      )
      stop(msg)
    }
    return(idx)
  }

  out <- c(0, 0, 0, 0)
  funs <- c("heter", "uncrt", "cat", "gp")
  types <- c("categorical", "categorical", "categorical", "continuous")
  fields <- c("", "", "cat_kernel", "gp_kernel")
  for (j in 1:4) {
    fun_covariate <- paste0(funs[j], "_covariate")
    names_type <- paste0(types[j], "_names")
    type_names <- nams[[names_type]]
    is_type <- !is.null(pf[[fun_covariate]])
    if (is_type) {
      cov_name <- pf[[fun_covariate]]
      fie <- fields[j]
      if (nchar(fie) == 0) {
        field <- funs[[j]]
      } else {
        field <- pf[[fie]]
      }
      idx <- check_type(cov_name, type_names, types[j], field)
    } else {
      idx <- 0
    }
    out[j] <- idx
  }

  # Return
  out
}

#' Helper for converting an lgpterm to numeric representation for Stan
#'
#' @param term an object of class \linkS4class{lgpterm}
#' @return a named list
check_term_factors <- function(term) {

  # Check for heter() expressions
  facs <- term@factors
  reduced <- reduce_factors_expr(facs, "heter")
  facs <- reduced$factors
  heter_covariate <- reduced$covariate

  # Check for uncrt() expressions
  reduced <- reduce_factors_expr(facs, "uncrt")
  facs <- reduced$factors
  uncrt_covariate <- reduced$covariate

  # Check compatibility
  if (!is.null(uncrt_covariate) && !is.null(heter_covariate)) {
    is_compatible <- (uncrt_covariate == heter_covariate)
    if (!is_compatible) {
      msg <- paste0(
        "Names of the covariates in uncrt() and heter() ",
        "expressions must match! Found = {",
        uncrt_covariate, ", ", heter_covariate, "}."
      )
      stop(msg)
    }
  }

  # Check for gp, gp_warp and gp_warp_vm expressions
  reduced <- reduce_factors_gp(facs)
  facs <- reduced$factors
  gp_covariate <- reduced$covariate
  gp_kernel <- reduced$kernel

  # Check for categorical kernel expression
  D <- length(facs)
  if (D == 0) {
    cat_covariate <- NULL
    cat_kernel <- NULL
  } else if (D == 1) {
    cat_covariate <- facs[[1]]@covariate
    cat_kernel <- facs[[1]]@fun
  } else {
    msg <- paste0(
      "Invalid term with expressions: ", as.character(term),
      ". Note that each term can contain at most one categ()",
      " or zerosum() expression."
    )
    stop(msg)
  }

  # Return a named list
  list(
    gp_covariate = gp_covariate,
    gp_kernel = gp_kernel,
    cat_covariate = cat_covariate,
    cat_kernel = cat_kernel,
    uncrt_covariate = uncrt_covariate,
    heter_covariate = heter_covariate
  )
}

#' Check for certain expressions in a term
#'
#' @param factors list of \linkS4class{lgpexpr} objects
#' @param expr the expression name to check
#' @return an updated list with no \code{expr} expressions, and name of
#' the covariate in the original \code{expr} expression
reduce_factors_expr <- function(factors, expr) {
  fun <- function(x) {
    x@fun == expr
  }
  idx <- which(sapply(factors, fun))
  H <- length(idx)
  if (H > 1) {
    msg <- paste0(
      "cannot have more than one '", expr, "' expression ",
      "in one term! found = ", H
    )
    stop(msg)
  }
  if (H == 1) {
    covariate <- factors[[idx]]@covariate
    factors[[idx]] <- NULL
  } else {
    covariate <- NULL
  }
  if (length(factors) < 1) {
    msg <- paste0(
      "there must be one <gp>, <gp_warp> or <gp_warp_vm> expression",
      " in each term involving the '", expr, "' expression!"
    )
    stop(msg)
  }
  list(factors = factors, covariate = covariate)
}

#' Check for gp expressions in a term
#'
#' @param factors list of \linkS4class{lgpexpr} objects
#' @return an updated list with no \code{gp*} expressions, and name of
#' the covariate in the original \code{gp*} expression
reduce_factors_gp <- function(factors) {
  fun <- function(x) {
    gp_names <- c("gp", "gp_warp", "gp_warp_vm")
    x@fun %in% gp_names
  }
  idx <- which(sapply(factors, fun))
  H <- length(idx)
  if (H > 1) {
    msg <- paste0(
      "cannot have more than one gp* expression ",
      "in one term! found = ", H
    )
    stop(msg)
  }
  if (H == 1) {
    covariate <- factors[[idx]]@covariate
    kernel <- factors[[idx]]@fun
    factors[[idx]] <- NULL
  } else {
    covariate <- NULL
    kernel <- NULL
  }
  list(factors = factors, covariate = covariate, kernel = kernel)
}
