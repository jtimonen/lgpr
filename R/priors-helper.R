#' Names of allowed  prior types
#'
#' @param idx an integer or NULL
#' @param return a character vector
prior_type_names <- function(idx = NULL) {
  names <- c(
    "Uniform", "Normal", "Student-t",
    "Gamma", "Inverse-Gamma", "Log-Normal"
  )
  names <- tolower(names)
  if (!is.null(idx)) {
    return(names[idx])
  } else {
    return(names)
  }
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
  names <- prior_type_names()
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
