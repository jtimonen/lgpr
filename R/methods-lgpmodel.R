# S4 METHODS --------------------------------------------------------------

#' @export
#' @describeIn lgpmodel Print information and summary about the object.
#' Returns \code{object} invisibly.
setMethod("show", "lgpmodel", function(object) {
  desc <- class_info("lgpmodel")
  cat(desc)
  cat("\n")
  model_summary(object)
})

#' @export
#' @describeIn lgpmodel Get a parameter summary (bounds and
#' priors). Returns a \code{data.frame}.
#' @param digits number of digits to show for floating point numbers
setMethod("parameter_info", "lgpmodel", function(object, digits = 3) {
  si <- get_stan_input(object)
  prior_to_df(si, digits = digits)
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
#' @describeIn lgpmodel Get names of model components.
setMethod("component_names", "lgpmodel", function(object) {
  rownames(get_component_encoding(object))
})

#' @export
#' @describeIn lgpmodel Determine if inference of the model requires sampling
#' the latent signal \code{f} (and its components).
setMethod("is_f_sampled", "lgpmodel", function(object) {
  object@sample_f
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

  # Helper function
  model_info <- function(object) {
    model <- object_to_model(object)
    stan_list <- get_stan_input(object)
    str1 <- as.character(model@model_formula)
    str2 <- likelihood_as_str(dollar(stan_list, "obs_model"))
    dat <- get_data(object)
    N <- nrow(dat)
    D <- ncol(dat)
    line1 <- paste0("Formula: ", str1)
    line2 <- paste0("Likelihood: ", str2)
    line3 <- paste0("Data: ", N, " observations, ", D, " variables")
    line4 <- approx_info(model)
    out <- paste0(line1, "\n", line2, "\n", line3, "\n", line4)
    return(out)
  }

  brief <- model_info(model)
  cat(brief)
  cat("Components:\n")
  print(component_info(model))
  cat("\n")
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
  print(parameter_info(model, digits))
  bi <- beta_teff_idx_info(model)
  beta_map <- dollar(bi, "beta")
  teff_map <- dollar(bi, "teff")
  het_z <- dollar(model@stan_input, "het_z")
  unc_z <- dollar(model@stan_input, "unc_z")
  if (!is.null(beta_map)) {
    cat(
      "\nConnection between categories of <",
      het_z, "> and beta parameter indices:\n",
      sep = ""
    )
    print(beta_map)
  }
  if (!is.null(teff_map)) {
    cat(
      "\nConnection between categories of <",
      unc_z, "> and teff parameter indices:\n",
      sep = ""
    )
    print(teff_map)
  }
  cat("\n")
  cat(misc_info(model), "\n")
  invisible(object)
}

#' @export
#' @rdname model_summary
param_summary <- function(object, digits = 3) {
  model <- object_to_model(object)
  parameter_info(model, digits)
}

# Info about mapping categories to parameter indices
beta_teff_idx_info <- function(object) {
  model <- object_to_model(object)
  beta_map <- dollar(model@stan_input, "BETA_IDX_MAP")
  teff_map <- dollar(model@stan_input, "TEFF_IDX_MAP")
  list(beta = beta_map, teff = teff_map)
}

# Information about used possible approximation
approx_info <- function(object) {
  model <- object_to_model(object)
  if (is_approximate(model)) {
    si <- get_stan_input(model)
    num_bf <- dollar(si, "num_bf")
    scale_bf <- dollar(si, "scale_bf")
    s1 <- paste0("num_bf = [", paste(num_bf, collapse = ", "), "]")
    s2 <- paste0("scale_bf = [", paste(scale_bf, collapse = ", "), "]")
    desc <- paste0("Approximation info: ", s1, ", ", s2)
  } else {
    desc <- ""
  }
  return(desc)
}

misc_info <- function(object) {
  model <- object_to_model(object)
  info <- model@info
  desc <- paste0(
    "Created on ", dollar(info, "created"), " with lgpr ",
    dollar(info, "lgpr_version"), "."
  )
  return(desc)
}

# Categorial covariate information
covariate_info.cat <- function(object) {
  model <- object_to_model(object)
  vn <- dollar(model@var_info, "var_names")
  nam <- dollar(vn, "z")
  if (is.null(nam)) {
    return(NULL)
  }
  num_levels <- dollar(model@stan_input, "Z_M")
  levels <- dollar(model@stan_input, "Z_levels")
  level_names <- c()
  J <- length(nam)
  for (j in seq_len(J)) {
    a <- levels[[j]]
    a <- if (length(a) > 4) "..." else paste(a, collapse = ", ")
    level_names[j] <- a
  }
  df <- data.frame(nam, num_levels, level_names)
  colnames(df) <- c("Factor", "#Levels", "Values")
  rownames(df) <- seq_len(dim(df)[1])
  return(df)
}

# Continuous covariate information
covariate_info.cont <- function(object) {
  model <- object_to_model(object)
  vn <- dollar(model@var_info, "var_names")
  nam <- dollar(vn, "x")
  if (is.null(nam)) {
    return(NULL)
  }
  mask <- dollar(model@stan_input, "X_mask")
  num_nan <- rowSums(mask == 1)
  df <- data.frame(nam, num_nan)
  colnames(df) <- c("Variable", "#Missing")
  rownames(df) <- seq_len(dim(df)[1])
  return(df)
}


# Determine if model is approximation
is_approximate <- function(object) {
  model <- object_to_model(object)
  si <- get_stan_input(model)
  num_bf <- dollar(si, "num_bf")
  if (any(num_bf > 0)) {
    return(TRUE)
  }
  FALSE
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

# Convert the Stan input encoding of a prior to a human-readable data frame
prior_to_df <- function(stan_input, digits = 3) {

  # Positive parameters
  check_positive(digits)
  pnames <- c("alpha", "ell", "wrp", "sigma", "phi")
  df <- NULL
  for (p in pnames) {
    df <- rbind(df, prior_to_df_pos(stan_input, p, digits))
  }

  # Beta
  num_het <- dollar(stan_input, "num_het")
  if (num_het > 0) {
    num_beta <- dollar(stan_input, "num_beta")
    df_bet <- prior_to_df_unit(stan_input, "beta", num_beta, digits)
    df <- rbind(df, df_bet)
  }

  # Gamma
  obs_model <- dollar(stan_input, "obs_model")
  if (obs_model == 5) {
    df_gam <- prior_to_df_unit(stan_input, "gamma", 1, digits)
    df <- rbind(df, df_gam)
  }

  # Effect time
  num_unc <- dollar(stan_input, "num_unc")
  if (num_unc > 0) {
    df_p <- prior_to_df_teff(stan_input, digits)
    df <- rbind(df, df_p)
  }

  return(df)
}

# Helper function for converting prior representation to human readable df
prior_to_df_pos <- function(stan_input, parname, digits) {
  prior <- dollar(stan_input, paste0("prior_", parname))
  hyper <- dollar(stan_input, paste0("hyper_", parname))
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
prior_to_df_unit <- function(stan_input, parname, num, digits) {
  hyper <- dollar(stan_input, paste0("hyper_", parname))
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
prior_to_df_teff <- function(stan_input, digits) {
  num_teff <- dollar(stan_input, "num_teff")
  prior <- dollar(stan_input, "prior_teff")
  type <- prior[1]
  backwards <- prior[2]
  hyper <- dollar(stan_input, "hyper_teff")
  zero <- dollar(stan_input, "teff_zero")
  lower <- dollar(stan_input, "teff_lb")
  upper <- dollar(stan_input, "teff_ub")
  pnames <- rep("foo", num_teff)
  dnames <- rep("foo", num_teff)
  bounds <- rep("foo", num_teff)
  for (j in seq_len(num_teff)) {
    par <- paste0("teff[", j, "]")
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

# Get the c_chat Stan input or a vector of zeros
get_chat <- function(object) {
  si <- get_stan_input(object)
  model <- object_to_model(object)
  if (is_f_sampled(model)) {
    c_hat <- dollar(si, "c_hat")
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
  if (is_f_sampled(object)) {
    stop(
      "Response variable is not normalized if f is sampled! Set ",
      "original = TRUE."
    )
  }
  si <- get_stan_input(object)
  out <- as.vector(dollar(si, "y"))
  return(out)
}

# Get response variable name
get_y_name <- function(object) {
  model <- object_to_model(object)
  dollar(model@var_names, "y")
}

# Get the Stan model used by a model
get_stan_model <- function(object) {
  model <- object_to_model(object)
  if (is_approximate(model)) {
    model_name <- "lgp_latent_bf"
  } else if (is_f_sampled(model)) {
    model_name <- "lgp_latent"
  } else {
    model_name <- "lgp"
  }
  stanmodels[[model_name]] # global variable (list of all pkg models)
}

# Get Stan input
get_stan_input <- function(object) {
  model <- object_to_model(object)
  return(model@stan_input)
}

# Get integer matrix encoding component types
get_component_encoding <- function(object) {
  si <- get_stan_input(object)
  dollar(si, "components")
}

# Get raw original data
get_data <- function(object) {
  model <- object_to_model(object)
  return(model@data)
}

# Get number of observations
get_num_obs <- function(object) {
  dollar(get_stan_input(object), "N")
}

# Get number of components
get_num_comps <- function(object) {
  dollar(get_stan_input(object), "J")
}

# Get observation model (human readable string)
get_obs_model <- function(object) {
  lh <- dollar(get_stan_input(object), "obs_model")
  likelihood_as_str(lh)
}

# Get number of trials (binomial or BB model)
get_num_trials <- function(object) {
  num_trials <- dollar(get_stan_input(object), "y_num_trials")
  as.vector(num_trials)
}
