# Create common Stan input needed for all models
create_model.base <- function(formula, data, options, prior, prior_only,
                              verbose) {

  # Data, formula and properties common for all models
  data <- convert_to_data_frame(data)
  lgp_formula <- parse_formula(formula, data, verbose)
  opts <- parse_options(options, prior_only, verbose)
  covariates <- parse_covariates(data, lgp_formula)
  components <- parse_components(lgp_formula, covariates)
  idx_maps <- create_index_maps(covariates, components)
  DIMS <- c(components, idx_maps)
  pri <- parse_prior(prior, DIMS, verbose)

  # Variable names
  var_names <- list(
    y = lgp_formula@y_name,
    x = rownames(dollar(covariates, "X")),
    z = rownames(dollar(covariates, "Z")),
    het_z = dollar(idx_maps, "het_z"),
    unc_z = dollar(idx_maps, "unc_z")
  )

  # Create the 'lgpmodel' object
  new("lgpmodel",
    model_formula = lgp_formula,
    raw_data = data,
    # Dimensions
    N = dollar(covariates, "N"),
    J = dollar(components, "J"),
    num_X = dollar(covariates, "num_X"),
    num_Z = dollar(covariates, "num_Z"),
    num_ell = dollar(components, "num_ell"),
    num_wrp = dollar(components, "num_wrp"),
    num_beta = dollar(idx_maps, "num_beta"),
    num_xpar = dollar(idx_maps, "num_xpar"),
    num_het = dollar(components, "num_het"),
    num_unc = dollar(components, "num_unc"),
    # Component encoding and param bounds
    components = dollar(components, "components"),
    # Covariates
    X = dollar(covariates, "X"),
    X_mask = dollar(covariates, "X_mask"),
    X_scale = dollar(covariates, "X_scale"),
    Z = dollar(covariates, "Z"),
    Z_M = dollar(covariates, "Z_M"),
    Z_levels = dollar(covariates, "Z_levels"),
    BETA_IDX = dollar(idx_maps, "BETA_IDX"),
    XPAR_IDX = dollar(idx_maps, "XPAR_IDX"),
    idx_maps = list(
      beta = dollar(idx_maps, "BETA_IDX_MAP"),
      xpar = dollar(idx_maps, "XPAR_IDX_MAP")
    ),
    # Options
    delta = dollar(opts, "delta"),
    vm_params = dollar(opts, "vm_params"),
    is_verbose = dollar(opts, "is_verbose"),
    is_likelihood_skipped = dollar(opts, "is_likelihood_skipped"),
    # Prior
    prior_alpha = dollar(pri, "prior_alpha"),
    prior_ell = dollar(pri, "prior_ell"),
    prior_wrp = dollar(pri, "prior_wrp"),
    prior_xpar = dollar(pri, "prior_xpar"),
    hyper_alpha = dollar(pri, "hyper_alpha"),
    hyper_ell = dollar(pri, "hyper_ell"),
    hyper_wrp = dollar(pri, "hyper_wrp"),
    hyper_xpar = dollar(pri, "hyper_xpar"),
    hyper_beta = dollar(pri, "hyper_beta"),
    xpar_zero = dollar(pri, "xpar_zero"),
    xpar_lb = dollar(pri, "xpar_lb"),
    xpar_ub = dollar(pri, "xpar_ub"),
    # Other
    var_names = var_names,
    info = creation_info()
  )
}

# Checks if formula is in advanced format and translates if not
parse_formula <- function(formula, data, verbose = FALSE) {
  advanced <- is_advanced_formula(formula)
  if (!advanced) formula <- formula_to_advanced(formula, data, verbose)
  fp <- as.character(formula)
  formula_str <- paste(fp[2], fp[1], fp[3])
  log_info(paste0("Formula interpreted as: ", formula_str), verbose)
  parse_formula_advanced(formula)
}

# Misc info for created objects
creation_info <- function() {
  list(
    created = date(),
    lgpr_version = utils::packageVersion("lgpr")
  )
}


# Parse given prior (common parameters)
parse_prior <- function(prior, DIMS, verbose) {
  num_unc <- dollar(DIMS, "num_unc")
  num_wrp <- dollar(DIMS, "num_wrp")
  par_names <- c("alpha", "ell", "wrp", "beta", "xpar", "xpar_info")
  filled <- fill_prior(prior, num_unc, par_names)
  defaulted <- defaulting_info(filled, verbose)
  wrp_defaulted <- "wrp" %in% defaulted
  if (num_wrp > 0 && wrp_defaulted) {
    model_desc <- "involves a gp_ns() or gp_vm() expression"
    msg <- warn_msg_default_prior("input warping steepness", "wrp", model_desc)
    warning(msg)
  }
  raw <- dollar(filled, "prior")
  parse_prior_common(raw, DIMS)
}


# Create mapping from observation index to index of beta or xpar parameter
create_index_maps <- function(covs, comps) {
  components <- dollar(comps, "components")
  Z <- dollar(covs, "Z")
  X_mask <- dollar(covs, "X_mask")
  Z_levs <- dollar(covs, "Z_levels")
  N <- ncol(Z)
  het <- create_expanding(components, Z, X_mask, "het()", 6, Z_levs)
  unc <- create_expanding(components, Z, X_mask, "unc()", 7, Z_levs)
  het <- reduce_index_maps(het, "het()", N, "beta")
  unc <- reduce_index_maps(unc, "unc()", N, "xpar")
  BETA_IDX <- dollar(het, "idx_expand")
  XPAR_IDX <- dollar(unc, "idx_expand")

  # Return
  list(
    BETA_IDX = BETA_IDX,
    XPAR_IDX = XPAR_IDX,
    BETA_IDX_MAP = dollar(het, "map"),
    XPAR_IDX_MAP = dollar(unc, "map"),
    num_beta = length(unique(BETA_IDX[BETA_IDX > 1])),
    num_xpar = length(unique(XPAR_IDX[XPAR_IDX > 1])),
    het_z = dollar(het, "covariate_name"),
    unc_z = dollar(unc, "covariate_name")
  )
}

# Parse the given common modeling options
parse_options <- function(options, prior_only, verbose) {

  # Default options
  input <- options
  opts <- list(delta = 1e-8, vm_params = c(0.025, 1))

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
  stan_switches <- list(
    is_likelihood_skipped = as.integer(prior_only),
    is_verbose = as.integer(verbose)
  )
  opts <- c(opts, stan_switches)
  return(opts)
}


# Create covariate data for Stan input
parse_covariates <- function(data, lgp_formula) {
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
    N = as.integer(N),
    num_X = as.integer(num_X),
    num_Z = as.integer(num_Z),
    Z = Z,
    Z_M = as.integer(Z_M),
    Z_levels = Z_levels,
    X = X,
    X_mask = X_mask,
    X_scale = X_scale
  )
}

# Create model components data for Stan input
parse_components <- function(model_formula, covariates) {
  components <- create_components_encoding(model_formula, covariates)
  list(
    J = nrow(components),
    components = components,
    num_ell = as.integer(sum(components[, 4] != 0)),
    num_wrp = as.integer(sum(components[, 5] != 0)),
    num_het = as.integer(sum(components[, 6] != 0)),
    num_unc = as.integer(sum(components[, 7] != 0))
  )
}
