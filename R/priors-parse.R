#' Parse given prior
#'
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters.
#' @param stan_input a list of stan input fields
#' @param obs_model observation model as integer
#' @return a named list of parsed options
parse_prior <- function(prior, stan_input, obs_model) {
  filled <- fill_prior(prior)
  str1 <- paste(filled$specified, collapse = ", ")
  str2 <- paste(filled$defaulted, collapse = ", ")
  msg1 <- paste0("User-specified priors found for: {", str1, "}")
  msg2 <- paste0("Default priors used for: {", str2, "}")
  info <- paste0(msg1, "\n", msg2, "\n")
  to_stan <- parse_prior_full(filled$prior, stan_input, obs_model)
  list(
    to_stan = to_stan,
    info = info
  )
}

#' Fill a partially defined prior
#'
#' @inheritParams parse_prior
#' @return a named list
fill_prior <- function(prior) {
  names <- c("alpha", "ell", "wrp", "sigma", "phi", 
             "beta", "effect_time")
  defaulted <- c()
  specified <- c()
  for (name in names) {
    if (name %in% names(prior)) {
      # User-specified prior
      specified <- c(specified, name)
      pr <- prior[[name]]
      is_named <- !is.null(names(pr))
      if (is_named) {
        pr <- list(pr)
      }
      prior[[name]] <- pr
    } else {
      # Default prior
      defaulted <- c(defaulted, name)
      prior[[name]] <- list(prior_default(name))
    }
  }
  list(
    prior = prior,
    defaulted = defaulted,
    specified = specified
  )
}


#' Parse given fully defined prior
#'
#' @inheritParams parse_prior
#' @return a named list of parsed options
parse_prior_full <- function(prior, stan_input, obs_model) {
  
  # Count number of different parameters
  nums <- stan_input[c("num_comps", "num_ell", "num_ns")]
  num_sigma <- as.numeric(obs_model == 1)
  num_phi <- as.numeric(obs_model == 3)
  num_uncrt <- stan_input$num_uncrt
  num_heter <- stan_input$num_heter
  num_cases <- stan_input$num_cases
  pnames <- c("alpha", "ell", "wrp", "sigma", "phi")
  nums <- c(unlist(nums), num_sigma, num_phi)
  names(nums)[4:5] <- c("num_sigma", "num_phi")
  
  # Common parameters
  common <- c()
  K <- length(pnames)
  for (k in seq_len(K)) {
    
    # Ensure that prior definition exists for this parameter
    desc <- field_must_exist(prior, pnames[k])
    
    # Prior type and hyper params to stan input format
    pp <- parse_prior_single(desc, nums[k])
    f1 <- paste0("prior_", pnames[k])
    f2 <- paste0("hyper_", pnames[k])
    common[[f1]] <- pp$prior
    common[[f2]] <- pp$hyper
  }
  
  # Beta parameters
  common[["hyper_beta"]] <- create_hyper_beta(prior, num_heter)

  # TODO: Change these
  other <- list(
    prior_teff = repvec(c(2, 0, 1), num_uncrt > 0),
    hyper_teff = repvec(c(0, 1, 0), num_uncrt > 0),
    
    teff_obs = repvec(rep(0, num_cases), num_uncrt > 0),
    teff_lb = repvec(rep(0, num_cases), num_uncrt > 0),
    teff_ub = repvec(rep(0, num_cases), num_uncrt > 0)
  )
  
  # Return
  c(common, other)
}

#' Create the hyper_beta input to Stan
#'
#' @inheritParams parse_prior_full
#' @param num_heter the \code{num_heter} input created for Stan
#' @return an array of size (0, 2) or (1, 2)
create_hyper_beta <- function(prior, num_heter){
  desc <- field_must_exist(prior, "beta")
  L <- length(desc)
  if (L != 1) {
    stop("length of prior$beta must be 1!")
  }
  bet <- c(desc[[1]]$alpha, desc[[1]]$beta)
  n_rows <- as.numeric(num_heter > 0)
  out <- repvec(bet, n_rows)
  return(out)
}

#' Parse the given prior for a single parameter type
#'
#' @param desc An unnamed list of prior descriptions.
#' @param num number of parameters of this type
#' @return a named list of parsed options
parse_prior_single <- function(desc, num) {
  nams <- names(desc)
  if (!is.null(nams)) {
    str <- paste(nams, collapse = ", ")
    stop("the list <desc> should not have names! found = {", str, "}")
  }
  L <- length(desc)
  if (L != num) {
    if (L == 1) {
      
      # The parameter type has the same prior in all components
      out <- prior_to_num(desc[[1]])
      prior <- repvec(out$prior, num)
      hyper <- repvec(out$hyper, num)
    } else {
      msg <- paste0(
        "<desc> should be a list of length 1 or ", num,
        "! found = ", L
      )
      stop(msg)
    }
  } else {
    
    # The parameter type has possibly different prior in different components
    prior <- repvec(c(0, 0), L)
    hyper <- repvec(c(0, 0, 0), L)
    for (j in seq_len(L)) {
      out <- prior_to_num(desc[[j]])
      prior[j, ] <- repvec(out$prior, num)
      hyper[j, ] <- repvec(out$hyper, num)
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
prior_to_num <- function(desc) {
  types <- prior_type_names()
  distribution_name <- desc$dist
  dist_num <- argument_check(distribution_name, types)
  fields <- names(desc)
  fields <- fields[!(fields %in% c("dist", "square"))]
  hyper <- position_hyper_params(desc[fields])
  list(
    prior = c(dist_num, as.numeric(desc$square)),
    hyper = hyper,
    hyper_names = fields
  )
}

#' Position the hyper parameters from a list to a vector that goes to Stan
#'
#' @param desc Hyperparameters as a named list
#' @return three real numbers
position_hyper_params <- function(desc) {
  hyper <- c(0, 0, 0)
  NAMES <- names(desc)
  for (name in NAMES) {
    val <- desc[[name]]
    H1 <- c("mu", "alpha", "nu")
    H2 <- c("sigma", "beta")
    if (name %in% H1) {
      hyper[1] <- val
    } else if (name %in% H2) {
      hyper[2] <- val
    } else {
      msg <- paste0(
        "invalid hyperparameter name '", name, "'! must be one of {",
        paste(c(H1, H2), collapse = ", ") , "}"
      )
      stop(msg)
    }
  }
  return(hyper)
}

