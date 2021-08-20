
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
