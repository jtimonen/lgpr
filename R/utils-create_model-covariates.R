#' Parse the covariates and model components from given data and formula
#'
#' @inheritParams create_model.likelihood
#' @param model_formula an object of class \linkS4class{lgpformula}
#' @param X_scale Scale factors for the continuous covariates
#' \itemize{
#'   \item an existing vector of scaling values
#'   \item \code{NA}, in which case such vector is create by computing
#'   the standard deviations from \code{data}
#' }
#' @return parsed input to Stan and covariate scaling, and other info
#' @family internal model creation functions
create_model.covs_and_comps <- function(data, model_formula,
                                        X_scale, verbose) {

  # Check that data is a data.frame and that all covariates exist in it
  NAMES <- unique(rhs_variables(model_formula@terms))
  log_progress("Parsing covariates and components...", verbose)
  check_df_with(data, NAMES)

  # Create the inputs to Stan
  covs <- stan_data_covariates(data, NAMES, X_scale)
  comps <- stan_data_components(model_formula, covs)
  print(comps)

  expanding <- stan_data_expanding(covs, comps)
  stan_data <- c(covs, comps, expanding)

  # Other info
  Z_levels <- dollar(covariates, "Z_levels")

  # Variable names
  var_names <- list(
    y = model_formula@y_name,
    x = rownames(dollar(covariates, "X")),
    z = rownames(dollar(covariates, "Z")),
  )

  # Return
  list(
    to_stan = to_stan,
    Z_levels = Z_levels,
    caseid_map = dollar(expanding, "caseid_map"),
    var_names = var_names
  )
}

# Create covariate data for Stan input
stan_data_covariates <- function(data, NAMES, X_scale) {
  check_unique(NAMES)
  check_not_null(X_scale)
  N <- dim(data)[1]
  scl_exists <- !is.na(X_scale)

  # Continuous
  X <- list()
  X_mask <- list()
  X_scale_new <- c()
  X_names <- c()
  num_X <- 0

  # Categorical
  Z <- list()
  Z_levels <- list()
  Z_M <- c()
  Z_names <- c()
  num_Z <- 0

  for (name in NAMES) {
    RAW <- data[[name]]
    if (is(RAW, "factor")) {

      # A categorical covariate
      num_Z <- num_Z + 1
      n_na <- sum(is.na(RAW))
      if (n_na > 0) {
        msg <- paste0(n_na, " missing values for factor '", name, "'!")
        stop(msg)
      }
      Z[[num_Z]] <- as.numeric(RAW)
      Z_M[num_Z] <- length(levels(RAW))
      Z_levels[[num_Z]] <- levels(RAW)
      Z_names[num_Z] <- name
    } else {

      # Continuous covariate, masking and scale
      num_X <- num_X + 1
      is_na <- is.na(RAW)
      X_mask[[num_X]] <- as.numeric(is_na) # store locations of NA, NaN
      X_NONAN <- RAW
      X_NONAN[is_na] <- 0 # create version where NA, NaN are replaced by 0
      X[[num_X]] <- X_NONAN
      X_names[num_X] <- name
      X_scale_new[num_X] <- stats::sd(RAW, na.rm = TRUE)
    }
  }

  # Convert to Stan input format
  Z <- list_to_matrix(Z, N)
  Z_M <- vector_to_array(Z_M)
  X <- list_to_matrix(X, N)
  X_mask <- list_to_matrix(X_mask, N)
  X_scale_new <- vector_to_array(X_scale_new)
  X_scale <- if (scl_exists) X_scale else X_scale_new

  # Name things
  names(Z_levels) <- Z_names
  names(Z_M) <- Z_names
  rownames(Z) <- Z_names
  rownames(X) <- X_names
  rownames(X_mask) <- X_names
  names(X_scale) <- X_names

  # Return list
  list(
    num_X = num_X,
    num_Z = num_Z,
    Z = Z,
    Z_M = Z_M,
    Z_levels = Z_levels,
    X = X,
    X_mask = X_mask,
    X_scale = X_scale
  )
}

# Helper function
vector_to_array <- function(x) {
  if (!is.null(x)) {
    x <- array(x, dim = length(x))
  } else {
    x <- array(0, dim = c(0))
  }
  return(x)
}

# Create model components data for Stan input
stan_data_components <- function(model_formula, covariates) {
  components <- stan_data_components.create(model_formula, covariates)
  list(
    J = nrow(components),
    components = components,
    num_ell = sum(components[, 4] != 0),
    num_wrp = sum(components[, 5] != 0),
    num_het = sum(components[, 6] != 0),
    idx_unc = 0 # TODO
  )
}

# Create an integer matrix that encodes component type info, and some other
# Stan inputs
stan_data_components.create <- function(model_formula, covariates) {
  terms <- model_formula@terms@summands
  J <- length(terms)
  comps <- array(0, dim = c(J, 7))
  for (j in seq_len(J)) {
    comps[j, ] <- term_to_numeric(terms[[j]], covariates)
  }
  colnames(comps) <- c("z", "k(z)", "x", "k(x)", "w", "h", "u")
  rownames(comps) <- term_names(model_formula@terms)
  as.matrix(comps)
}

# Create mapping from observation index to index of beta or teff parameter
stan_data_expanding <- function(covs, comps) {
  components <- dollar(comps, "components")
  Z <- dollar(covs, "Z")
  X_mask <- dollar(covs, "X_mask")
  lst <- stan_data_expanding.create(components, Z, X_mask)
  idx_expand <- dollar(lst, "idx_expand")
  num_bt <- length(unique(idx_expand[idx_expand > 1]))

  # Return
  list(
    idx_expand = idx_expand,
    num_bt = num_bt,
    caseid_map = dollar(lst, "map")
  )
}

# Map a list of terms to their "names"
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

# An lgpterm to numeric representation (6 integers) for Stan
term_to_numeric <- function(term, covariates) {
  opts <- rep(0, 7)
  parsed <- parse_term_factors(term)
  Z_names <- rownames(dollar(covariates, "Z"))
  X_names <- rownames(dollar(covariates, "X"))
  opts[1:2] <- term_to_numeric.z(parsed, Z_names)
  opts[3:5] <- term_to_numeric.x(parsed, X_names)
  opts[6] <- term_to_numeric.het(parsed, Z_names)
  opts[7] <- term_to_numeric.unc(parsed, Z_names)
  return(opts)
}

# Create options [1-2] to 'components' Stan input
term_to_numeric.z <- function(parsed, Z_names) {
  opts <- c(0, 0)
  ker_name <- dollar(parsed, "z_kernel")
  if (!is.null(ker_name)) {

    # Covariate
    cn <- dollar(parsed, "z_covariate")
    type <- "categorical"
    opts[1] <- term_to_numeric.check_cov_type(cn, Z_names, type, ker_name)

    # Kernel
    opts[2] <- check_allowed(ker_name, allowed = c("categ", "zs"))
  }
  return(opts)
}

# Create options [3-5] to 'components' Stan input
term_to_numeric.x <- function(parsed, X_names) {
  ker_name <- dollar(parsed, "x_kernel")
  opts <- c(0, 0, 0)
  if (!is.null(ker_name)) {

    # Covariate
    cn <- dollar(parsed, "x_covariate")
    type <- "continuous"
    opts[1] <- term_to_numeric.check_cov_type(cn, X_names, type, ker_name)

    # Kernel
    opts[2] <- 1 # all three use EQ kernel
    idx_k <- check_allowed(ker_name, allowed = c("gp", "gp_ns", "gp_vm"))
    opts[3] <- idx_k - 1
  }
  return(opts)
}

# Create option [6] to 'components' Stan input
term_to_numeric.het <- function(parsed, Z_names) {
  opt <- 0
  cov_name <- dollar(parsed, "het_covariate")
  if (!is.null(cov_name)) {
    type <- "categorical"
    opt <- term_to_numeric.check_cov_type(cov_name, Z_names, type, "het")
  }
  return(opt)
}

# Create option [7] to 'components' Stan input
term_to_numeric.unc <- function(parsed, Z_names) {
  opt <- 0
  cov_name <- dollar(parsed, "unc_covariate")
  if (!is.null(cov_name)) {
    type <- "categorical"
    opt <- term_to_numeric.check_cov_type(cov_name, Z_names, type, "unc")
  }
  return(opt)
}

# Helper function
term_to_numeric.check_cov_type <- function(cov_name, allowed, type, fun) {
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

# Helper for converting an lgpterm to numeric representation for Stan
parse_term_factors <- function(term) {

  # Remove all het() expressions and store covariate name
  facs <- term@factors
  reduced <- reduce_factors_expr(facs, "het")
  facs <- dollar(reduced, "factors")
  het_covariate <- dollar(reduced, "covariate")

  # Remove all unc() expressions and store covariate name
  reduced <- reduce_factors_expr(facs, "unc")
  facs <- dollar(reduced, "factors")
  unc_covariate <- dollar(reduced, "covariate")

  # Remove for gp(), gp_ns(), or gp_vm() expression and store covariate name
  # and kernel name
  reduced <- reduce_factors_gp(facs)
  facs <- dollar(reduced, "factors")
  x_covariate <- dollar(reduced, "covariate")
  x_kernel <- dollar(reduced, "kernel")

  # Remove for zs() or categ() expression and store covariate name
  # and kernel name
  D <- length(facs)
  if (D == 0) {
    z_covariate <- NULL
    z_kernel <- NULL
  } else if (D == 1) {
    z_covariate <- facs[[1]]@covariate
    z_kernel <- facs[[1]]@fun
  } else {
    msg <- paste0(
      "Invalid term with expressions: ", as.character(term),
      ". \nNote that each term can contain\n",
      " * at most one het() or unc() expression,\n",
      " * at most one gp(), gp_ns() or gp_vm() expression AND\n",
      " * at most one zs() or categ() expression.\n",
    )
    stop(msg)
  }

  # Return a named list
  list(
    x_covariate = x_covariate,
    x_kernel = x_kernel,
    z_covariate = z_covariate,
    z_kernel = z_kernel,
    unc_covariate = unc_covariate,
    het_covariate = het_covariate
  )
}

# Check for certain expressions in a formula term
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
      "there must be one gp(), gp_ns() or gp_vm() expression",
      " in each term involving the '", expr, "' expression!"
    )
    stop(msg)
  }
  list(factors = factors, covariate = covariate)
}

# Check for gp, gp_ns and gp_vm expressions in a term
reduce_factors_gp <- function(factors) {
  fun <- function(x) {
    gp_names <- c("gp", "gp_ns", "gp_vm")
    x@fun %in% gp_names
  }
  idx <- which(sapply(factors, fun))
  H <- length(idx)
  if (H > 1) {
    msg <- paste0(
      "Cannot have more than one gp(), gp_ns() or gp_vm() expression ",
      "in one term (found = ", H,
      ")! If you used the | syntax in your model formula, make sure that",
      " the covariate on the left side of | is continuous and the one on",
      " the right side is categorical."
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

# Creating the idx_expand input for Stan
stan_data_expanding.create <- function(components, x_cat, x_cont_mask) {
  pick <- create_idx_expand_picker(components, x_cat)
  x_fac <- dollar(pick, "x_fac")
  factor_name <- dollar(pick, "factor_name")
  inds <- which(components[, 4] + components[, 7] > 0)
  inds <- as.numeric(inds) # remove names
  L <- length(inds)
  if (L == 0) {
    idx_expand <- x_fac
    map <- NULL
  } else {
    i_cont <- components[inds, 9]
    to_reduce <- x_cont_mask[i_cont, ]
    if (length(i_cont) == 1) {
      to_reduce <- repvec(to_reduce, 1)
    }
    idx_mask <- reduce_rows(to_reduce)
    map <- map_factor_to_caseid(x_fac, idx_mask, factor_name)
    idx_expand <- map_caseid_to_row(x_fac, map) + 1 # note the plus 1
  }

  # Return
  list(map = map, idx_expand = idx_expand)
}

# Helper for creating the idx_expand input for Stan
create_idx_expand_picker <- function(components, x_cat) {
  n_obs <- dim(x_cat)[2]
  inds <- c(components[, 4], components[, 7])
  inds <- as.numeric(inds[inds != 0])
  J <- length(inds)
  if (J == 0) {
    x_fac <- rep(1, n_obs)
    name <- NULL
  } else {
    all_same <- all(inds == inds[1])
    if (!all_same) {
      str <- paste(inds, collapse = ", ")
      msg <- paste0(
        "The het() and unc() expressions must have the same ",
        "categorical covariate in every term! ",
        "Found inds = {", str, "}"
      )
      stop(msg)
    }
    i1 <- inds[1]
    x_fac <- as.numeric(x_cat[i1, ])
    name <- rownames(x_cat)[i1]
  }

  # Return
  list(x_fac = x_fac, factor_name = name)
}


# Helper for creating the idx_expand input for Stan
map_factor_to_caseid <- function(x_fac, x_cont_mask, factor_name) {
  id <- c()
  case_id <- c()
  num_cases <- 0
  facs <- unique(x_fac)
  for (u in facs) {
    inds <- which(x_fac == u)
    vals <- x_cont_mask[inds]
    L <- length(vals)
    all0 <- all.equal(vals, rep(0, L)) == TRUE
    all1 <- all.equal(vals, rep(1, L)) == TRUE
    if (all0) {
      # case individual, set its case_id
      num_cases <- num_cases + 1
      id[num_cases] <- u
      case_id[num_cases] <- num_cases
    } else if (all1) {
      # do nothing, not a case
    } else {
      msg <- paste0(
        "inconsistent x_cont_mask values for observations where ",
        factor_name, " = ", u
      )
      stop(msg)
    }
  }
  out <- data.frame(id, case_id)
  colnames(out) <- c(factor_name, "case_id")
  return(out)
}

# Helper for creating the idx_expand input for Stan
map_caseid_to_row <- function(x_fac, map) {
  out <- rep(0, length(x_fac))
  num_levels <- dim(map)[1]
  for (j in seq_len(num_levels)) {
    id <- map[j, 1]
    caseid <- map[j, 2]
    inds <- which(x_fac == id)
    out[inds] <- caseid
  }
  return(out)
}
