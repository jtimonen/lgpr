# TERM TO NUMERIC ---------------------------------------------------------

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


# PARSE TERM FACTORS ------------------------------------------------------

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
      " * at most one het() expression,\n",
      " * at most one unc() expression,\n",
      " * at most one gp(), gp_ns() or gp_vm() expression AND\n",
      " * at most one zs() or categ() expression.\n"
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
