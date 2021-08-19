# MAIN LEVEL --------------------------------------------------------------

#' Parse the covariates and model components from given data and formula
#'
#' @inheritParams create_model.likelihood
#' @param model_formula an object of class \linkS4class{lgpformula}
#' @return parsed input to Stan and covariate scaling, and other info
#' @family internal model creation functions
create_model.covs_and_comps <- function(data, model_formula, verbose) {

  # Check that data is a data.frame and that all covariates exist in it
  NAMES <- unique(rhs_variables(model_formula@terms))
  log_progress("Parsing covariates and components...", verbose)
  check_df_with(data, NAMES)

  # Create the inputs to Stan
  covs <- stan_data_covariates(data, NAMES)
  comps <- stan_data_components(model_formula, covs)
  expanding <- stan_data_expanding(covs, comps)
  stan_data <- c(covs, comps, expanding)

  # Variable names
  var_names <- list(
    y = model_formula@y_name,
    x = rownames(dollar(covs, "X")),
    z = rownames(dollar(covs, "Z"))
  )

  # Return
  list(to_stan = stan_data, var_names = var_names)
}


# COVARIATES --------------------------------------------------------------

# Create covariate data for Stan input
stan_data_covariates <- function(data, NAMES) {
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
stan_data_components <- function(model_formula, covariates) {
  components <- stan_data_components.create(model_formula, covariates)
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
stan_data_components.create <- function(model_formula, covariates) {
  terms <- model_formula@terms@summands
  J <- length(terms)
  comps <- array(0, dim = c(J, 7))
  for (j in seq_len(J)) {
    comps[j, ] <- term_to_numeric(terms[[j]], covariates)
  }
  colnames(comps) <- c("iz", "k(z)", "ix", "k(x)", "wrp", "het", "unc")
  rownames(comps) <- stan_data_components.term_names(model_formula@terms)
  as.matrix(comps)
}

# Map a list of terms to their "names"
stan_data_components.term_names <- function(rhs) {
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
stan_data_expanding <- function(covs, comps) {
  components <- dollar(comps, "components")
  Z <- dollar(covs, "Z")
  X_mask <- dollar(covs, "X_mask")
  Z_levs <- dollar(covs, "Z_levels")
  N <- ncol(Z)
  het <- stan_data_expanding.create(components, Z, X_mask, "het()", 6, Z_levs)
  unc <- stan_data_expanding.create(components, Z, X_mask, "unc()", 7, Z_levs)
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
stan_data_expanding.create <- function(components, Z, X_mask, expr, icol,
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
    maps <- stan_data_expanding.maps(components, j, icol, X_mask, Z, Z_levels)
    lst[[cntr]] <- maps
    names(lst)[cntr] <- rownames(components)[j]
  }
  list(lst = lst, covariate_name = unames)
}

# Index maps for component j
stan_data_expanding.maps <- function(components, j, icol, X_mask, Z, Z_levels) {
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
  is_missing <- which(x_mask == 1)
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
