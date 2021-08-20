#' Parse given prior
#'
#' @inheritParams create_model.likelihood
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters. See the "Defining priors" section below
#' (\code{\link{lgp}}).
#' @param stan_input a list of 'Stan' input data fields
#' @return a named list of parsed options
#' @family internal model creation functions
create_model.prior <- function(prior, stan_input, verbose) {
  log_progress("Parsing prior...", verbose)
  num_unc <- dollar(stan_input, "num_unc")
  num_wrp <- dollar(stan_input, "num_wrp")
  par_names <- c("alpha", "ell", "wrp", "sigma", "phi", "beta", "gamma")
  par_names <- c(par_names, "effect_time", "effect_time_info")
  filled <- fill_prior(prior, num_unc)
  spec <- dollar(filled, "specified")
  dflt <- dollar(filled, "defaulted")
  str1 <- paste(spec, collapse = ", ")
  str2 <- paste(dflt, collapse = ", ")
  wrp_defaulted <- "wrp" %in% dflt
  if (num_wrp > 0 && wrp_defaulted) {
    model_desc <- "involves a gp_ns() or gp_vm() expression"
    msg <- warn_msg_default_prior("input warping steepness", "wrp", model_desc)
    warning(msg)
  }
  msg1 <- paste0("User-specified priors found for: {", str1, "}.")
  msg2 <- paste0(
    "If any of the following parameters are included in the",
    " model, default priors are used for them: {", str2, "}."
  )
  info <- paste0(msg1, "\n", msg2, "\n")
  log_info(info, verbose)

  raw <- dollar(filled, "prior")
  to_stan <- parse_prior_full(raw, stan_input)
  list(
    to_stan = to_stan,
    raw = raw
  )
}
