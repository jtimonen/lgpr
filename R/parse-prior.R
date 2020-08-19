#' Parse the given prior
#'
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters. It is recommended to first create this using the function
#' \code{\link{prior_default}}, and then possibly modify it.
#' @param stan_input a list of stan input fields
#' @param obs_model observation model as integer
#' @return a named list of parsed options
parse_prior <- function(prior, stan_input, obs_model) {
  num_comps <- stan_input$num_comps
  num_ell <- stan_input$num_ell
  num_ns <- stan_input$num_ns
  num_heter <- stan_input$num_heter
  num_uncrt <- stan_input$num_uncrt
  num_cases <- stan_input$num_cases
  num_sigma <- as.numeric(obs_model == 1)
  num_phi <- as.numeric(obs_model == 3)

  list(
    prior_alpha = repvec(c(3, 0), num_comps),
    prior_ell = repvec(c(6, 0), num_ell),
    prior_wrp = repvec(c(6, 0), num_ns),
    prior_sigma = repvec(c(6, 0), num_sigma),
    prior_phi = repvec(c(6, 0), num_phi),
    prior_teff = repvec(c(2, 0, 1), num_uncrt > 0),

    hyper_alpha = repvec(c(20, 0, 0), num_comps),
    hyper_ell = repvec(c(1, 1, 0), num_ell),
    hyper_wrp = repvec(c(1, 1, 0), num_ns),
    hyper_sigma = repvec(c(1, 1, 0), num_sigma),
    hyper_phi = repvec(c(1, 1, 0), num_phi),
    hyper_beta = repvec(c(0.2, 0.2), num_heter > 0),
    hyper_teff = repvec(c(0, 1, 0), num_uncrt > 0),

    teff_obs = repvec(rep(0, num_cases), num_uncrt > 0),
    teff_lb = repvec(rep(0, num_cases), num_uncrt > 0),
    teff_ub = repvec(rep(0, num_cases), num_uncrt > 0)
  )
}

#' Convert the given prior in a Stan input to a human-readable format
#'
#' @inheritParams prior_to_df_oneparam
#' @return a data frame
prior_to_df <- function(stan_input, digits = 3) {
  pnames <- c("alpha", "ell", "wrp", "sigma", "phi")
  df <- NULL
  for (p in pnames) {
    df_p <- prior_to_df_oneparam(stan_input, p, digits)
    if (is.null(df)) {
      df <- df_p
    } else {
      df <- rbind(df, df_p)
    }
  }
  return(df)
}

#' Convert the given prior for one parameter type into a human-readable format
#'
#' @param stan_input a list of stan input fields
#' @inheritParams prior_to_char
#' @return a data frame
prior_to_df_oneparam <- function(stan_input, parname, digits) {
  prior <- stan_input[[paste0("prior_", parname)]]
  hyper <- stan_input[[paste0("hyper_", parname)]]
  D <- dim(prior)[1]
  pnames <- rep("foo", D)
  dnames <- rep("foo", D)
  for (j in seq_len(D)) {
    par <- paste0(parname, "[", j, "]")
    out <- prior_to_char(par, prior[j, ], hyper[j, ], digits)
    pnames[j] <- out$parname
    dnames[j] <- out$distribution
  }
  df <- data.frame(pnames, dnames)
  colnames(df) <- c("Parameter", "Prior")
  return(df)
}

#' Human-readable prior statement
#'
#' @param parname parameter name
#' @param prior two integers
#' @param hyper three real numbers
#' @param digits number of digits to show for floats
#' @return A list.
prior_to_char <- function(parname, prior, hyper, digits) {
  hyper <- round(hyper, digits)

  # Check distribution type
  tp <- prior[1]
  if (tp < 1 || tp > 6) {
    stop("Prior type must 1, 2, 3, 4, 5 or 6! Found = ", tp)
  }
  names <- c(
    "Uniform", "Normal", "Student-t",
    "Gamma", "Inverse-Gamma", "Log-Normal"
  )
  pname <- names[tp]

  # Check if there is a transform
  tf <- prior[2]
  if (tf == 1) {
    parname <- paste0("(", parname, ")^2")
  } else if (tf == 0) {
    parname <- parname
  } else {
    msg <- paste0(
      "Transform must be either 0 (identity) or 1 (squaring). ",
      "Found = ", tf
    )
    stop(msg)
  }

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
