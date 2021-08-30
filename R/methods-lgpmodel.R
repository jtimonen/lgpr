# S4 METHODS --------------------------------------------------------------

#' @export
#' @describeIn lgpmodel Print information and summary about the object.
#' Returns \code{object} invisibly.
setMethod("show", "lgpmodel", function(object) {
  desc <- class_info(class(object))
  cat(desc)
  cat("\n")
  model_summary(object)
})

#' @export
#' @describeIn lgpmodel Get a parameter summary (bounds and
#' priors). Returns a \code{data.frame}.
#' @param digits number of digits to show for floating point numbers
setMethod("parameter_info", "lgpmodel", function(object, digits = 3) {
  check_positive(digits)
  prior_to_df.common(object, digits)
})

#' @export
#' @describeIn MarginalGPModel Get a parameter summary (bounds and
#' priors). Returns a \code{data.frame}.
#' @param digits number of digits to show for floating point numbers
setMethod("parameter_info", "MarginalGPModel", function(object, digits = 3) {
  check_positive(digits)
  df_common <- prior_to_df.common(object, digits)
  df_add <- prior_to_df.marginal(object, digits)
  rbind(df_common, df_add)
})

#' @export
#' @describeIn LatentGPModel Get a parameter summary (bounds and
#' priors). Returns a \code{data.frame}.
#' @param digits number of digits to show for floating point numbers
setMethod("parameter_info", "LatentGPModel", function(object, digits = 3) {
  check_positive(digits)
  df_common <- prior_to_df.common(object, digits)
  df_add <- prior_to_df.latent(object, digits)
  rbind(df_common, df_add)
})

#' @export
#' @describeIn lgpmodel Get a data frame with information about each model
#' component.
setMethod("component_info", "lgpmodel", function(object) {
  comps <- get_component_encoding(object)
  nam <- colnames(comps)
  df <- data.frame(comps)
  colnames(df) <- nam
  return(df)
})

#' @export
#' @describeIn lgpmodel Get covariate information.
setMethod("covariate_info", "lgpmodel", function(object) {
  info1 <- covariate_info.cont(object)
  info2 <- covariate_info.cat(object)
  list(continuous = info1, categorical = info2)
})

#' @export
#' @describeIn lgpmodel Get number of model components.
setMethod("num_components", "lgpmodel", function(object) {
  object@J
})

#' @export
#' @describeIn lgpmodel Get names of model components.
setMethod("component_names", "lgpmodel", function(object) {
  rownames(get_component_encoding(object))
})

#' @export
#' @describeIn lgpmodel Determine if inference of the model requires sampling
#' the latent signal \code{f} (and its components).
setMethod("is_f_sampled", "lgpmodel", function(object) {
  !is_f_marginalized(object)
})

# Is f marginalized?
is_f_marginalized <- function(object) {
  m <- object_to_model(object)
  if (isa(m, "MarginalGPModel")) {
    return(TRUE)
  }
  FALSE
}

# Is the response variable normalized during inference?
is_y_normalized <- function(object) {
  om <- get_obs_model(object)
  return(om == "gaussian")
}

# Is the model approximate?
is_approx <- function(object) {
  m <- object_to_model(object)
  if (isa(m, "LatentGPModelApprox")) {
    return(TRUE)
  }
  FALSE
}

#' @export
#' @describeIn MarginalGPModel Get name of corresponding 'Stan' model.
setMethod("get_stanmodel", "MarginalGPModel", function(object) {
  stanmodels[["lgp_marginal"]]
})

#' @export
#' @describeIn LatentGPModel Get name of corresponding 'Stan' model.
setMethod("get_stanmodel", "LatentGPModel", function(object) {
  stanmodels[["lgp_latent"]]
})

#' @export
#' @describeIn LatentGPModelApprox Get name of corresponding 'Stan' model.
setMethod("get_stanmodel", "LatentGPModelApprox", function(object) {
  stanmodels[["lgp_latent_approx"]]
})


# SUMMARY AND INFO --------------------------------------------------------

#' Print a model summary.
#'
#' @export
#' @param object a model or fit
#' @param digits number of digits to round floats to
#' @return \code{object} invisibly.
model_summary <- function(object, digits = 3) {
  model <- object_to_model(object)
  misc_summary(model)
  model_summary_brief(model)
  component_summary(model)
  covariate_summary(model)
  parameter_summary(model, digits)
  invisible(object)
}

# Print brief model summary
model_summary_brief <- function(object) {
  model <- object_to_model(object)
  str <- as.character(model@model_formula)
  N <- get_num_obs(model)
  cat("Data:", N, "observations\n")
  cat("Normalize response:", is_y_normalized(object), "\n")
  cat("Marginalize f:", is_f_marginalized(object), "\n")
  cat("Approximation:", is_approx(object), "\n")
  cat("Formula:", str, "\n")
}

# Print component summary
component_summary <- function(model) {
  cat("Components (", num_components(model), "):\n", sep = "")
  print(component_info(model))
  cat("\n")
}

# Print covariate summary
covariate_summary <- function(model) {
  ci <- covariate_info(model)
  info_cont <- dollar(ci, "continuous")
  info_cat <- dollar(ci, "categorical")
  if (!is.null(info_cont)) {
    print(info_cont)
    cat("\n")
  }
  if (!is.null(info_cat)) {
    print(info_cat)
    cat("\n")
  }
}

# Print a parameter summary of a model
parameter_summary <- function(model, digits = 3) {
  print(parameter_info(model, digits))
  bi <- beta_xpar_idx_info(model)
  beta_map <- dollar(bi, "beta")
  xpar_map <- dollar(bi, "xpar")
  if (!is.null(beta_map)) {
    het_z <- dollar(model@var_names, "het_z")
    cat(
      "\nConnection between categories of <",
      het_z, "> and beta parameter indices:\n",
      sep = ""
    )
    print(beta_map)
  }
  if (!is.null(xpar_map)) {
    unc_z <- dollar(model@var_names, "unc_z")
    cat(
      "\nConnection between categories of <",
      unc_z, "> and xpar parameter indices:\n",
      sep = ""
    )
    print(xpar_map)
  }
  cat("\n")
}

# Print information about used possible approximation
approx_summary <- function(object) {
  model <- object_to_model(object)
  num_bf <- model@num_bf
  scale_bf <- model@scale_bf
  s1 <- paste0("num_bf = [", paste(num_bf, collapse = ", "), "]")
  s2 <- paste0("scale_bf = [", paste(scale_bf, collapse = ", "), "]")
  desc <- paste0("Approximation info: ", s1, ", ", s2, "\n")
  return(desc)
}

# Print misc info about a model
misc_summary <- function(object) {
  model <- object_to_model(object)
  info <- model@info
  desc <- paste0("Created with lgpr ", dollar(info, "lgpr_version"), ".\n")
  cat(desc)
}


# Categorial covariate information
covariate_info.cat <- function(object) {
  model <- object_to_model(object)
  vn <- model@var_names
  nam <- dollar(vn, "z")
  if (is.null(nam)) {
    return(NULL)
  }
  num_levels <- model@Z_M
  levels <- model@Z_levels
  level_names <- c()
  J <- length(nam)
  for (j in seq_len(J)) {
    a <- levels[[j]]
    a <- if (length(a) > 4) "..." else paste(a, collapse = ", ")
    level_names[j] <- a
  }
  df <- data.frame(num_levels, level_names)
  colnames(df) <- c("#Categories", "Values")
  rownames(df) <- nam
  return(df)
}

# Continuous covariate information
covariate_info.cont <- function(object) {
  model <- object_to_model(object)
  vn <- model@var_names
  nam <- dollar(vn, "x")
  if (is.null(nam)) {
    return(NULL)
  }
  mask <- model@X_mask
  num_nan <- rowSums(mask == 1)
  df <- data.frame(num_nan)
  colnames(df) <- c("#Missing")
  rownames(df) <- nam
  return(df)
}

# Info about mapping categories to parameter indices
beta_xpar_idx_info <- function(object) {
  object <- object_to_model(object)
  beta_map <- dollar(object@idx_maps, "beta")
  xpar_map <- dollar(object@idx_maps, "xpar")
  list(beta = beta_map, xpar = xpar_map)
}

# Helper function for plots
#
# @param object model or fit
# @param x x-axis variable name
# @param group_by grouping variable name (use \code{NULL} for no grouping)
# @return a data frame
create_plot_df <- function(object, x = "age", group_by = "id") {

  # Get x-axis variable
  dat <- get_data(object)
  x_name <- x
  x <- dollar(dat, x_name)
  check_type(x, "numeric")

  # Get grouping factor
  x_grp <- create_grouping_factor(dat, group_by) # util

  # Get response
  y <- get_y(object, original = TRUE)
  y_name <- get_y_name(object)
  df <- data.frame(x_grp, x, y)
  group_by <- if (is.na(group_by)) "group__" else group_by
  colnames(df) <- c(group_by, x_name, y_name)
  return(df)
}


# PRIOR TO DF -------------------------------------------------------------

# Convert parsed prior to a human-readable data frame
prior_to_df.common <- function(model, digits = 3) {

  # Positive parameters
  check_positive(digits)
  pnames <- c("alpha", "ell", "wrp")
  df <- NULL
  for (p in pnames) {
    df <- rbind(df, prior_to_df_pos(model, p, digits))
  }

  # Beta
  num_het <- model@num_het
  if (num_het > 0) {
    num_beta <- model@num_beta
    df_bet <- prior_to_df_unit(model, "beta", num_beta, digits)
    df <- rbind(df, df_bet)
  }

  # Effect time
  num_unc <- model@num_unc
  if (num_unc > 0) {
    df_p <- prior_to_df_xpar(model, digits)
    df <- rbind(df, df_p)
  }

  return(df)
}

# Convert parsed prior to a human-readable data frame
prior_to_df.marginal <- function(model, digits = 3) {
  df <- prior_to_df_pos(model, "sigma", digits)
  return(df)
}

# Convert parsed prior to a human-readable data frame
prior_to_df.latent <- function(model, digits = 3) {
  LH <- model@obs_model
  if (LH == 1) {
    df <- prior_to_df_pos(model, "sigma", digits)
  } else if (LH == 3) {
    df <- prior_to_df_pos(model, "phi", digits)
  } else if (LH == 5) {
    df <- prior_to_df_unit(model, "gamma", 1, digits)
  } else {
    df <- NULL
  }
  return(df)
}

# Helper function for converting prior representation to human readable df
prior_to_df_pos <- function(model, parname, digits) {
  prior <- slot(model, paste0("prior_", parname))
  hyper <- slot(model, paste0("hyper_", parname))
  D <- dim(prior)[1]
  pnames <- rep("foo", D)
  dnames <- rep("foo", D)
  bounds <- rep("foo", D)
  for (j in seq_len(D)) {
    par <- paste0(parname, "[", j, "]")
    out <- prior_to_str(par, prior[j, ], hyper[j, ], digits)
    tpar <- dollar(out, "parname")
    pnames[j] <- par
    dnames[j] <- paste0(tpar, " ~ ", dollar(out, "distribution"))
    bounds[j] <- "[0, Inf)"
  }
  df <- data.frame(pnames, bounds, dnames)
  colnames(df) <- c("Parameter", "Bounds", "Prior")
  return(df)
}

# Helper function for converting prior representation to human readable df
prior_to_df_unit <- function(model, parname, num, digits) {
  hyper <- slot(model, paste0("hyper_", parname))
  check_positive(digits)
  a <- round(hyper[1], digits = digits)
  b <- round(hyper[2], digits = digits)
  dist <- paste0("beta(", a, ", ", b, ")")
  bounds <- "[0, 1]"
  nam1 <- paste0(parname, "[1]")
  nam2 <- paste0(parname, "[1-", num, "]")
  par <- if (num > 1) nam2 else nam1
  dist <- paste0(par, " ~ ", dist)
  df <- data.frame(par, bounds, dist)
  colnames(df) <- c("Parameter", "Bounds", "Prior")
  return(df)
}

# Helper function for converting prior representation to human readable df
prior_to_df_xpar <- function(model, digits) {
  num_xpar <- model@num_xpar
  if (num_xpar == 0) {
    return(NULL)
  }
  prior <- model@prior_xpar
  type <- prior[1]
  backwards <- prior[2]
  hyper <- model@hyper_xpar
  zero <- model@xpar_zero
  lower <- model@xpar_lb
  upper <- model@xpar_ub
  pnames <- rep("foo", num_xpar)
  dnames <- rep("foo", num_xpar)
  bounds <- rep("foo", num_xpar)
  for (j in seq_len(num_xpar)) {
    par <- paste0("xpar[", j, "]")
    tpar <- par
    tpar <- minus.append(tpar, zero[j])
    tpar <- minus.prepend(tpar, backwards)
    out <- prior_to_str(par, c(type, 0), hyper, digits)
    pnames[j] <- par
    dnames[j] <- paste0(tpar, " ~ ", dollar(out, "distribution"))
    bounds[j] <- paste0("[", lower[j], ", ", upper[j], "]")
  }
  df <- data.frame(pnames, bounds, dnames)
  colnames(df) <- c("Parameter", "Bounds", "Prior")
  return(df)
}

# Append minus and val to a string if val is not zero
minus.append <- function(str, val) {
  if (val != 0) str <- paste0(str, " - ", val)
  return(str)
}

# Prepend minus to a string
minus.prepend <- function(str, prepend) {
  if (prepend != 0) str <- paste0(" - (", str, ")")
  return(str)
}

# Human-readable prior statement
prior_to_str <- function(parname, prior, hyper, digits) {
  hyper <- round(hyper, digits)

  # Check distribution type
  tp <- prior[1]
  check_allowed(tp, seq_len(6))
  names <- prior_type_names()
  pname <- names[tp]

  # Check if there is a transform
  tf <- prior[2]
  check_allowed(tf, c(0, 1))
  if (tf == 1) parname <- paste0("(", parname, ")^2")

  # Get prior statement
  if (tp %in% c(2, 4, 5, 6)) {
    str <- paste0(pname, "(", hyper[1], ",", hyper[2], ")")
  } else if (tp == 3) {
    str <- paste0(pname, "(", hyper[1], ")")
  } else {
    str <- paste(pname, sep = "")
  }

  # Return
  list(parname = parname, distribution = str)
}


# GETTERS -----------------------------------------------------------------

# Get the c_hat vector a vector of zeros
get_chat <- function(object) {
  model <- object_to_model(object)
  if (is_f_sampled(model)) {
    c_hat <- as.vector(model@c_hat)
  } else {
    N <- get_num_obs(model)
    c_hat <- rep(0.0, N)
  }
  return(c_hat)
}

# Get response variable measurements on original or normalized scale
get_y <- function(object, original = TRUE) {
  if (original) {
    y_name <- get_y_name(object)
    dat <- get_data(object)
    return(dollar(dat, y_name))
  }
  if (!is_y_normalized(object)) {
    stop(
      "Response variable is not normalized! Set ",
      "original = TRUE."
    )
  }
  out <- as.vector(object_to_model(object)@y)
  return(out)
}

# Get response variable name
get_y_name <- function(object) {
  model <- object_to_model(object)
  dollar(model@var_names, "y")
}

# Get integer matrix encoding component types
get_component_encoding <- function(object) {
  object_to_model(object)@components
}

# Get raw original data
get_data <- function(object) {
  object_to_model(object)@raw_data
}

# Get number of observations
get_num_obs <- function(object) {
  object_to_model(object)@N
}

# Get observation model (human readable string)
get_obs_model <- function(object) {
  model <- object_to_model(object)
  if (isa(model, "MarginalGPModel")) {
    LH <- 1L
  } else {
    LH <- model@obs_model
  }
  likelihood_as_str(LH)
}

# Get number of trials (binomial or BB model)
get_num_trials <- function(object) {
  as.vector(object_to_model(object)@y_num_trials)
}

# CREATING STAN DATA ------------------------------------------------------

# Create Stan data list for a model
create_standata <- function(model) {
  slots <- standata_slotnames(model)
  dat <- list()
  j <- 0
  for (s in slots) {
    j <- j + 1
    dat[[j]] <- slot(model, s)
  }
  names(dat) <- slots
  return(dat)
}

# Get names of needed slots for Stan data list
standata_slotnames <- function(model) {

  # Names defined in inst/stan/_common/data-general.stan
  common_general <- c(
    "is_verbose", "is_likelihood_skipped", "N", "J", "num_X", "num_Z",
    "num_ell", "num_wrp", "num_beta", "num_xpar", "num_het", "num_unc",
    "components", "delta", "vm_params"
  )

  # Names defined in inst/stan/_common/data-covariates.stan
  common_covariates <- c(
    "X", "X_mask", "X_scale", "Z", "Z_M", "BETA_IDX", "XPAR_IDX"
  )

  # Names defined in inst/stan/_common/data-priors.stan
  common_priors <- c(
    "prior_alpha", "prior_ell", "prior_wrp", "prior_xpar",
    "hyper_alpha", "hyper_ell", "hyper_wrp", "hyper_beta",
    "hyper_xpar", "xpar_zero", "xpar_lb", "xpar_ub"
  )

  slots <- c(common_general, common_covariates, common_priors)
  if (isa(model, "MarginalGPModel")) {
    # Names defined in inst/stan/lgp_marginal.stan
    slots_add <- c("y", "prior_sigma", "hyper_sigma")
  } else if (isa(model, "LatentGPModel")) {
    # Names defined in inst/stan/_latent/data-*.stan
    slots_add <- c(
      "y_int", "y", "c_hat", "y_num_trials", "obs_model",
      "prior_sigma", "prior_phi", "hyper_sigma", "hyper_phi", "hyper_gamma"
    )
  } else {
    slots_add <- c("y_int", "y", "c_hat", "y_num_trials")
  }
  slots <- c(slots, slots_add)
  return(slots)
}
