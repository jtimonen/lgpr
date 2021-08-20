# FORMULA -----------------------------------------------------------------

# Check if formula is in advanced format
is_advanced_formula <- function(formula) {
  f <- split_formula(formula)
  rhs <- dollar(f, "right")
  out <- grepl("(", rhs, fixed = TRUE)
  return(out)
}

# Split formula into left and right-hand side
split_formula <- function(formula) {
  check_type(formula, "formula")
  f_str <- as.character(formula)
  check_length(f_str, 3)
  out <- list(left = f_str[2], right = f_str[3])
  return(out)
}

# Translate formula to advanced syntax format
formula_to_advanced <- function(formula, data, verbose) {
  f <- split_formula(formula)
  rhs <- dollar(f, "right")
  y_name <- dollar(f, "left")
  rhs <- rhs_to_advanced(rhs, data, verbose, y_name)
  f_str <- paste(y_name, "~", rhs)
  out <- stats::formula(f_str)
  return(out)
}

# Translate formula right-hand side
rhs_to_advanced <- function(rhs, data, verbose, y_name) {
  types <- data_types(data, y_name, verbose)
  terms <- rhs_to_terms(rhs)
  rhs_adv <- ""
  idx <- 0
  for (t in terms) {
    idx <- idx + 1
    covs <- strsplit(t, split = "|", fixed = TRUE)[[1]]
    if (length(covs) == 1) {
      check_in_data(covs[1], data, "data")
      t1 <- dollar(types, covs[1])
      if (t1 == "numeric") {
        f1 <- "gp"
      } else if (t1 == "factor") {
        f1 <- "zs"
      } else {
        stop("Invalid data type at this point. Please report a bug.")
      }
      term <- enclose_fun(covs[1], f1)
    } else {
      check_length(covs, 2)
      check_in_data(covs[1], data, "data")
      check_in_data(covs[2], data, "data")
      e1 <- enclose_fun(covs[1], "gp")
      t2 <- dollar(types, covs[2])
      f2 <- if ("factor" %in% t2) "zs" else "gp"
      e2 <- enclose_fun(covs[2], f2)
      term <- paste0(e1, "*", e2)
    }
    rhs_adv <- if (idx == 1) term else paste(rhs_adv, "+", term)
  }
  return(rhs_adv)
}

# Split a formula right-hand side to terms
rhs_to_terms <- function(rhs) {
  x <- simplify_str(rhs)
  out <- strsplit(x, split = "+", fixed = TRUE)[[1]]
  return(out)
}

# Parse a model formula that uses the advanced formula syntax and
# return an object of class \linkS4class{lgpformula}
parse_formula_advanced <- function(formula) {
  f <- split_formula(formula)
  rhs <- dollar(f, "right")
  terms <- parse_rhs(rhs)
  y_name <- dollar(f, "left")
  new("lgpformula",
    y_name = y_name,
    terms = terms,
    call = paste(y_name, "~", rhs)
  )
}

# Parse a string representation of the right-hand side of an advanced formula
parse_rhs <- function(rhs) {
  terms <- rhs_to_terms(rhs)
  D <- length(terms)
  f <- parse_term(terms[1])
  if (D == 1) {
    out <- new("lgprhs", summands = list(f))
    return(out)
  } else {
    for (j in 2:D) {
      f <- f + parse_term(terms[j])
    }
    return(f)
  }
}

# Parse a string representation of one advanced formula term
parse_term <- function(term) {
  factors <- strsplit(term, split = "*", fixed = TRUE)[[1]] # split to factors
  D <- length(factors)
  if (D > 4) {
    stop("One formula term can have at most four expressions! found = ", D)
  }
  e <- parse_expr(factors[[1]])
  f <- new("lgpterm", factors = list(e))
  if (D == 1) {
    return(f)
  } else {
    for (j in 2:D) {
      e <- parse_expr(factors[[j]])
      f <- f * new("lgpterm", factors = list(e))
    }
    return(f)
  }
}

# Parse a string representation of one advanced formula expression
parse_expr <- function(expr) {
  num_open <- lengths(regmatches(expr, gregexpr("[(]", expr)))
  num_close <- lengths(regmatches(expr, gregexpr("[)]", expr)))
  if (num_open == 0) {
    msg <- paste0(
      "Zero opening brackets found in the expression <",
      expr, ">. Be sure not to mix the advanced and simple",
      " formula syntaxes. See ?lgp."
    )
    stop(msg)
  }
  ok_paren <- (num_open == 1) && (num_close == 1)
  if (!ok_paren) {
    msg <- paste0(
      "Each expression must contain exactly one opening and ",
      "closing parenthesis! Found ", num_open, " opening and ",
      num_close, " closing brackets in the expression <", expr, ">."
    )
    stop(msg)
  }
  parsed <- strsplit(expr, "[()]")[[1]]
  check_length(parsed, 2)
  out <- new("lgpexpr", covariate = parsed[2], fun = parsed[1])
  return(out)
}

# Get names of all variables appearing in a term of advanced formula
term_variables <- function(term) {
  a <- character()
  for (f in term@factors) {
    a <- c(a, f@covariate)
  }
  return(a)
}

# Get names of all variables appearing on advanced formula right-hand side
rhs_variables <- function(rhs) {
  a <- character()
  for (s in rhs@summands) {
    a <- c(a, term_variables(s))
  }
  return(a)
}


# COMMON OPTIONS ------------------------------------------------------------

# Parse the given common modeling options
standata_common_options <- function(options, prior_only) {

  # Default options
  input <- options
  opts <- list(
    delta = 1e-8,
    vm_params = c(0.025, 1)
  )

  # Replace defaults if found from input
  for (opt_name in names(opts)) {
    if (opt_name %in% names(input)) {
      opts[[opt_name]] <- input[[opt_name]]
    }
  }

  # Validate and format for Stan input
  delta <- dollar(opts, "delta")
  vm_params <- dollar(opts, "vm_params")
  check_positive(delta)
  check_length(vm_params, 2)
  check_positive_all(vm_params)
  check_all_leq(vm_params, c(1, 1))

  # Return full options
  opts$is_likelihood_skipped <- as.numeric(prior_only)
  return(opts)
}


# COVARIATES --------------------------------------------------------------

# Create covariate data for Stan input
standata_covariates <- function(data, lgp_formula) {
  NAMES <- unique(rhs_variables(lgp_formula@terms))
  check_df_with(data, NAMES)
  check_unique(NAMES)
  N <- dim(data)[1]

  # Continuous
  X <- list()
  X_mask <- list()
  X_scale <- c()
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
      x_sd <- stats::sd(RAW, na.rm = TRUE)
      if (x_sd == 0) {
        stop("the variable <", name, "> has zero variance")
      }
      X_scale[num_X] <- x_sd
    }
  }

  # Convert to Stan input format
  Z <- list_to_matrix(Z, N)
  Z_M <- vector_to_array(Z_M)
  X <- list_to_matrix(X, N)
  X_mask <- list_to_matrix(X_mask, N)
  X_scale <- vector_to_array(X_scale)

  # Name things
  names(Z_levels) <- Z_names
  names(Z_M) <- Z_names
  rownames(Z) <- Z_names
  rownames(X) <- X_names
  rownames(X_mask) <- X_names
  names(X_scale) <- X_names

  # Return list
  list(
    N = N,
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


# COMPONENTS --------------------------------------------------------------

# Create model components data for Stan input
standata_components <- function(model_formula, covariates) {
  components <- create_components_encoding(model_formula, covariates)
  list(
    J = nrow(components),
    components = components,
    num_ell = sum(components[, 4] != 0),
    num_wrp = sum(components[, 5] != 0),
    num_het = sum(components[, 6] != 0),
    num_unc = sum(components[, 7] != 0)
  )
}

# Create an integer matrix that encodes component type info, and some other
# Stan inputs
create_components_encoding <- function(model_formula, covariates) {
  terms <- model_formula@terms@summands
  J <- length(terms)
  comps <- array(0, dim = c(J, 7))
  for (j in seq_len(J)) {
    comps[j, ] <- term_to_numeric(terms[[j]], covariates)
  }
  colnames(comps) <- c("iz", "k(z)", "ix", "k(x)", "wrp", "het", "unc")
  rownames(comps) <- term_names(model_formula@terms)
  as.matrix(comps)
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


# EXPANDING ---------------------------------------------------------------

# Create mapping from observation index to index of beta or teff parameter
standata_expanding <- function(covs, comps) {
  components <- dollar(comps, "components")
  Z <- dollar(covs, "Z")
  X_mask <- dollar(covs, "X_mask")
  Z_levs <- dollar(covs, "Z_levels")
  N <- ncol(Z)
  het <- standata_expanding.create(components, Z, X_mask, "het()", 6, Z_levs)
  unc <- standata_expanding.create(components, Z, X_mask, "unc()", 7, Z_levs)
  het <- reduce_index_maps(het, "het()", N, "beta")
  unc <- reduce_index_maps(unc, "unc()", N, "teff")
  BETA_IDX <- dollar(het, "idx_expand")
  TEFF_IDX <- dollar(unc, "idx_expand")

  # Return
  list(
    BETA_IDX = BETA_IDX,
    TEFF_IDX = TEFF_IDX,
    BETA_IDX_MAP = dollar(het, "map"),
    TEFF_IDX_MAP = dollar(unc, "map"),
    num_beta = length(unique(BETA_IDX[BETA_IDX > 1])),
    num_teff = length(unique(TEFF_IDX[TEFF_IDX > 1])),
    het_z = dollar(het, "covariate_name"),
    unc_z = dollar(unc, "covariate_name")
  )
}

# Helper function
standata_expanding.create <- function(components, Z, X_mask, expr, icol,
                                      Z_levels) {
  inds <- which(components[, icol] > 0) # indices of needed components
  inside_hets <- components[inds, icol]
  het_names <- rownames(Z)[inside_hets]
  unames <- unique(het_names)
  if (length(unames) > 1) {
    msg <- paste0(
      "Covariate inside the ", expr, " expression must be the same",
      " in every term that has such expression. Found  = {",
      paste(het_names, collapse = ", "), "}."
    )
    stop(msg)
  }
  lst <- list()
  cntr <- 0
  for (j in inds) {
    cntr <- cntr + 1
    maps <- standata_expanding.maps(components, j, icol, X_mask, Z, Z_levels)
    lst[[cntr]] <- maps
    names(lst)[cntr] <- rownames(components)[j]
  }
  list(lst = lst, covariate_name = unames)
}

# Index maps for component j
standata_expanding.maps <- function(components, j, icol, X_mask, Z, Z_levels) {
  Z_names <- rownames(Z)
  idx_x <- components[j, 3]
  idx_z <- components[j, icol]
  z <- Z[idx_z, ]
  x_mask <- X_mask[idx_x, ]
  map <- category_inds_to_param_inds(z, x_mask)
  colnames(map) <- dollar(Z_levels, Z_names[idx_z])
  idx_expand <- obs_inds_to_param_inds(z, map)
  list(
    idx_expand = idx_expand,
    map = map
  )
}

# Create a map from category to parameter index
category_inds_to_param_inds <- function(z, x_mask) {
  uval <- sort(unique(z))
  Q <- length(uval)
  param_idx <- rep(0, Q)
  pidx <- 0
  for (q in seq_len(Q)) {
    inds <- which(z == uval[q])
    if (any(x_mask[inds] == 0)) {
      pidx <- pidx + 1
      param_idx[q] <- pidx
    }
  }
  arr <- rbind(uval, param_idx)
  rownames(arr) <- c("category_idx", "param_idx")
  return(arr)
}

# Create a map from category to parameter index
obs_inds_to_param_inds <- function(z, map) {
  N <- length(z)
  idx_expand <- rep(1, N)
  for (n in seq_len(N)) {
    idx_expand[n] <- map[2, z[n]]
  }
  return(idx_expand + 1) # plus 1, because 1 means no param
}

# Check that the param_idx maps are same for each term, and reduce maps to one
reduce_index_maps <- function(idx_maps, expr, N, param_name) {
  maps <- dollar(idx_maps, "lst")
  cn <- dollar(idx_maps, "covariate_name")
  if (length(maps) == 0) {
    idx_expand <- array(0, dim = c(0, N)) # empty version
    map <- NULL
  } else {
    a <- sapply(maps, "[[", "idx_expand") # a is array of shape c(N, J)
    u <- unique(a, MARGIN = 2) # take unique columns
    if (ncol(u) > 1) {
      msg <- paste0(
        "If there are multiple terms that contain the ", expr,
        " expression, their continuous covariates have to either",
        " all be missing or all be available for any given data row."
      )
      stop(msg)
    }
    idx_expand <- array(as.vector(u), dim = c(1, N))
    map <- dollar(maps[[1]], "map")
    rownames(map)[1] <- paste0(cn, "_idx")
    rownames(map)[2] <- paste0(param_name, "_idx")
  }

  # Return
  list(
    idx_expand = idx_expand,
    map = map,
    covariate_name = cn
  )
}


# PRIOR -------------------------------------------------------------------

# Parse given prior (common parameters)
standata_common_prior <- function(prior, stan_input, verbose) {
  num_unc <- dollar(stan_input, "num_unc")
  num_wrp <- dollar(stan_input, "num_wrp")
  par_names <- c(
    "alpha", "ell", "wrp", "beta",
    "effect_time", "effect_time_info"
  )
  filled <- fill_prior(prior, num_unc, par_names)
  defaulted <- defaulting_info(filled, verbose)
  wrp_defaulted <- "wrp" %in% defaulted
  if (num_wrp > 0 && wrp_defaulted) {
    model_desc <- "involves a gp_ns() or gp_vm() expression"
    msg <- warn_msg_default_prior("input warping steepness", "wrp", model_desc)
    warning(msg)
  }
  raw <- dollar(filled, "prior")
  parse_prior_common(raw, stan_input)
}

# Names of allowed prior types
prior_type_names <- function(idx = NULL) {
  names <- c(
    "Uniform",
    "Normal",
    "Student-t",
    "Gamma",
    "Inv-Gamma",
    "Log-Normal"
  )
  names <- tolower(names)
  out <- if (!is.null(idx)) names[idx] else names
  return(out)
}

# Fill a partially defined prior
fill_prior <- function(prior, num_uncrt, par_names) {
  defaulted <- c()
  specified <- c()
  for (name in par_names) {
    if (name %in% names(prior)) {
      # User-specified prior
      specified <- c(specified, name)
      pr <- dollar(prior, name)
      prior[[name]] <- list_if_named(pr)
    } else {
      # Default prior
      defaulted <- c(defaulted, name)
      prior[[name]] <- default_prior(name, num_uncrt)
    }
  }
  list(
    prior = prior,
    defaulted = defaulted,
    specified = specified
  )
}

# Information about which parameters have default priors
defaulting_info <- function(filled, verbose) {
  spec <- dollar(filled, "specified")
  dflt <- dollar(filled, "defaulted")
  str1 <- paste(spec, collapse = ", ")
  str2 <- paste(dflt, collapse = ", ")
  msg1 <- paste0("User-specified priors found for: {", str1, "}.")
  msg2 <- paste0(
    "If any of the following parameters are included in the",
    " model, default priors are used for them: {", str2, "}."
  )
  info <- paste0(msg1, "\n", msg2, "\n")
  log_info(info, verbose)
  return(dflt)
}

# Parse common parameters of a prior
parse_prior_common <- function(prior, si) {
  DIM_BETA <- as.numeric(dollar(si, "num_het") > 0)
  DIM_TEFF <- as.numeric(dollar(si, "num_unc") > 0)
  num_teff <- dollar(si, "num_teff")
  alpha <- parse_prior_pos(prior, dollar(si, "J"), "alpha")
  ell <- parse_prior_pos(prior, dollar(si, "num_ell"), "ell")
  wrp <- parse_prior_pos(prior, dollar(si, "num_wrp"), "wrp")
  beta <- parse_prior_unit(prior, DIM_BETA, "beta")
  teff <- parse_prior_teff(prior, DIM_TEFF, num_teff)
  c(alpha, ell, wrp, beta, teff)
}

# Parse prior of a positive parameter
parse_prior_pos <- function(prior, DIM, par_name) {
  desc <- dollar(prior, par_name) # no [[1]]
  parsed <- list()
  pp <- parse_prior_single(desc, DIM)
  f1 <- paste0("prior_", par_name)
  f2 <- paste0("hyper_", par_name)
  parsed[[f1]] <- dollar(pp, "prior")
  parsed[[f2]] <- dollar(pp, "hyper")
  return(parsed)
}

# Parse a unit interval parameter prior (beta prior)
parse_prior_unit <- function(prior, DIM, par_name) {
  parsed <- list()
  desc <- dollar(prior, par_name)[[1]]
  hyper <- c(dollar(desc, "alpha"), dollar(desc, "beta"))
  f <- paste0("hyper_", par_name)
  parsed[[f]] <- repvec(hyper, DIM)
  return(parsed)
}

# Parse effect time parameter prior
parse_prior_teff <- function(prior, DIM, num_teff) {
  effect_time_info <- dollar(prior, "effect_time_info")[[1]]
  is_backwards <- as.numeric(dollar(effect_time_info, "backwards"))
  lower <- dollar(effect_time_info, "lower")
  upper <- dollar(effect_time_info, "upper")
  zero <- dollar(effect_time_info, "zero")
  lower <- ensure_len(lower, num_teff)
  upper <- ensure_len(upper, num_teff)
  zero <- ensure_len(zero, num_teff)

  desc <- dollar(prior, "effect_time")[[1]]
  out <- prior_to_num(desc)
  type <- dollar(out, "prior")
  hyper <- dollar(out, "hyper")
  prior <- c(type[1], is_backwards)

  # Return
  list(
    prior_teff = repvec(prior, DIM),
    hyper_teff = repvec(hyper, DIM),
    teff_zero = repvec(zero, DIM),
    teff_lb = repvec(lower, DIM),
    teff_ub = repvec(upper, DIM)
  )
}

# Parse given prior for a single parameter type
parse_prior_single <- function(desc, num) {
  check_not_named(desc)
  L <- length(desc)
  err_msg <- paste0("<desc> should have length 1 or ", num, "! found = ", L)
  if (L != num) {
    if (L == 1) {
      # The parameter type has the same prior in all components
      out <- prior_to_num(desc[[1]])
      prior <- repvec(dollar(out, "prior"), num)
      hyper <- repvec(dollar(out, "hyper"), num)
    } else {
      stop(err_msg)
    }
  } else {

    # The parameter type has possibly different prior in different components
    prior <- repvec(c(0, 0), L)
    hyper <- repvec(c(0, 0, 0), L)
    for (j in seq_len(L)) {
      desc_j <- desc[[j]]
      out <- prior_to_num(desc_j)
      pr <- repvec(dollar(out, "prior"), 1)
      hp <- repvec(dollar(out, "hyper"), 1)
      prior[j, ] <- pr
      hyper[j, ] <- hp
    }
  }
  list(
    prior = prior,
    hyper = hyper
  )
}

#' Convert given prior to numeric format
#'
#' @param desc Prior description as a named list, containing fields
#' \itemize{
#'   \item \code{dist} - Distribution name. Must be one of
#'   {'uniform', 'normal', 'student-t', 'gamma', 'inv-gamma', 'log-normal'}
#'   (case-insensitive)
#'   \item \code{square} - Is the prior for a square-transformed parameter.
#' }
#' Other list fields are interpreted as hyperparameters.
#' @return a named list of parsed options
prior_to_num <- function(desc) {
  types <- prior_type_names()
  distribution_name <- dollar(desc, "dist")
  dist_num <- check_allowed(distribution_name, types)
  fields <- names(desc)
  fields <- fields[!(fields %in% c("dist", "square"))]
  hyper <- position_hyper_params(desc[fields])
  sq <- dollar(desc, "square")
  list(
    prior = c(dist_num, as.numeric(sq)),
    hyper = hyper,
    hyper_names = fields
  )
}

# Position the hyper parameters from a list to a vector that goes to Stan
position_hyper_params <- function(desc) {
  hyper <- c(0, 0, 0)
  NAMES <- names(desc)
  H1 <- c("mu", "alpha", "nu")
  H2 <- c("sigma", "beta")
  for (name in NAMES) {
    check_allowed(name, c(H1, H2))
    val <- dollar(desc, name)
    idx <- if (name %in% H2) 2 else 1
    hyper[idx] <- val
  }
  return(hyper)
}

# Default prior, given parameter name and number of uncertain components
default_prior <- function(name, num_uncrt = NULL) {
  if (name == "effect_time_info") {
    desc <- default_prior_effect_time_info(num_uncrt)
  } else {
    desc <- list(default_priors_most(name))
  }
  return(desc)
}

# Default priors for most parameters
default_priors_most <- function(name) {
  allowed <- c(
    "alpha", # 1
    "ell", # 2
    "sigma", # 3
    "phi", # 4
    "wrp", # 5
    "beta", # 6
    "gamma", # 7
    "effect_time" # 8
  )
  idx <- check_allowed(name, allowed)
  if (idx == 1) {
    prior <- student_t(nu = 20) # alpha
  } else if (idx == 2) {
    prior <- log_normal(mu = 0, sigma = 1) # ell
  } else if (idx == 3) {
    prior <- igam(shape = 2, scale = 1, square = TRUE) # sigma
  } else if (idx == 4) {
    prior <- log_normal(mu = 1, sigma = 1, square = TRUE) # phi
  } else if (idx == 5) {
    prior <- igam(shape = 14, scale = 5) # wrp
  } else if (idx == 6) {
    prior <- bet(a = 0.2, b = 0.2) # beta
  } else if (idx == 7) {
    prior <- bet(a = 1, b = 1) # gamma
  } else {
    prior <- uniform() # effect_time
  }
  return(prior)
}

# Default effect time info
default_prior_effect_time_info <- function(num_uncrt) {
  check_not_null(num_uncrt)
  if (num_uncrt > 0) {
    # There is no default prior for effect time
    # TODO: more informative message?
    stop("you must specify 'effect_time_info' in <prior>!")
  } else {
    # Will not be used
    desc <- list(backwards = FALSE, lower = NaN, upper = NaN, zero = NaN)
    desc <- list(desc)
  }
  return(desc)
}
