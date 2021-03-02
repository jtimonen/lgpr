#' Parse a fully defined prior
#'
#' @inheritParams parse_prior
#' @name parse_prior_full
#' @return a named list to be given to Stan
#' @family prior parsers
NULL

#' @name parse_prior_full
parse_prior_full <- function(prior, stan_input, obs_model) {
  list_pos <- parse_prior_full_pos(prior, stan_input, obs_model)
  list_unit <- parse_prior_full_unit(prior, stan_input, obs_model)
  list_teff <- parse_prior_full_teff(prior, stan_input)
  c(list_pos, list_unit, list_teff)
}

#' @rdname parse_prior_full
parse_prior_full_pos <- function(prior, stan_input, obs_model) {
  par_names <- c("alpha", "ell", "wrp", "sigma", "phi")
  nums <- stan_input[c("num_comps", "num_ell", "num_ns")]
  num_sigma <- as.numeric(obs_model == 1)
  num_phi <- as.numeric(obs_model == 3)
  nums <- c(unlist(nums), num_sigma, num_phi)

  parsed <- c()
  K <- length(par_names)
  for (k in seq_len(K)) {
    desc <- dollar(prior, par_names[k])
    pp <- parse_prior_single(desc, nums[k])
    f1 <- paste0("prior_", par_names[k])
    f2 <- paste0("hyper_", par_names[k])
    parsed[[f1]] <- dollar(pp, "prior")
    parsed[[f2]] <- dollar(pp, "hyper")
  }
  return(parsed)
}


#' @rdname parse_prior_full
parse_prior_full_unit <- function(prior, stan_input, obs_model) {
  par_names <- c("beta", "gamma")
  num_gamma <- as.numeric(obs_model == 5)
  num_heter <- dollar(stan_input, "num_heter")
  num_heter <- as.numeric(num_heter > 0)
  nums <- c(num_heter, num_gamma)

  parsed <- c()
  K <- length(par_names)
  for (k in seq_len(K)) {
    desc <- dollar(prior, par_names[k])[[1]]
    hyper <- c(dollar(desc, "alpha"), dollar(desc, "beta"))
    f <- paste0("hyper_", par_names[k])
    parsed[[f]] <- repvec(hyper, nums[k])
  }
  return(parsed)
}

#' @rdname parse_prior_full
parse_prior_full_teff <- function(prior, stan_input) {
  num_bt <- dollar(stan_input, "num_bt")
  num_uncrt <- dollar(stan_input, "num_uncrt")

  effect_time_info <- dollar(prior, "effect_time_info")[[1]]
  is_backwards <- as.numeric(dollar(effect_time_info, "backwards"))
  lower <- dollar(effect_time_info, "lower")
  upper <- dollar(effect_time_info, "upper")
  zero <- dollar(effect_time_info, "zero")
  lower <- ensure_len(lower, num_bt)
  upper <- ensure_len(upper, num_bt)
  zero <- ensure_len(zero, num_bt)
  DIM <- as.numeric(num_uncrt > 0)

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

#' Parse the given prior for a single parameter type
#'
#' @param desc An unnamed list of prior descriptions.
#' @param num number of parameters of this type
#' @return a named list of parsed options
#' @family prior parsers
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

#' Parse the given prior to numeric format
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
#' @family prior parsers
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

#' Position the hyper parameters from a list to a vector that goes to Stan
#'
#' @param desc Hyperparameters as a named list
#' @return three real numbers
#' @family prior parsers
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
