#' Parse given prior
#'
#' @param prior A named list, defining the prior distribution of model
#' (hyper)parameters.
#' @param stan_input a list of stan input fields
#' @param obs_model observation model as integer
#' @return a named list of parsed options
parse_prior <- function(prior, stan_input, obs_model) {
  prior <- fill_prior(prior)
  parse_prior_full(prior, stan_input, obs_model)
}

#' Fill a partially defined prior
#'
#' @inheritParams parse_prior
#' @return a named list
fill_prior <- function(prior) {
  if (!is.null(prior)) {
    stop("Prior was not NULL")
  }
  st20 <- student_t(nu = 20)
  ln11 <- log_normal(mu = 1, sigma = 1)
  changeThis <- normal(mu = 1, sigma = 0.5)
  list(
    alpha = list(st20),
    ell = list(ln11),
    wrp = list(changeThis),
    phi = list(ln11),
    sigma = list(ln11)
  )
}

#' Parse given fully defined prior
#'
#' @inheritParams parse_prior
#' @return a named list of parsed options
parse_prior_full <- function(prior, stan_input, obs_model) {
  nums <- stan_input[c("num_comps", "num_ell", "num_ns")]
  num_sigma <- as.numeric(obs_model == 1)
  num_phi <- as.numeric(obs_model == 3)
  num_heter <- stan_input$num_heter
  num_uncrt <- stan_input$num_uncrt
  num_cases <- stan_input$num_cases

  pnames <- c("alpha", "ell", "wrp", "sigma", "phi")
  nums <- c(unlist(nums), num_sigma, num_phi)
  names(nums)[4:5] <- c("num_sigma", "num_phi")

  common <- c()
  K <- length(pnames)
  for (k in seq_len(K)) {
    name <- pnames[k]
    desc <- prior[[name]]
    if (is.null(desc)) {
      stop("<prior> must contain a field named ", name)
    }
    pp <- parse_prior_single(desc, nums[k])
    f1 <- paste0("prior_", name)
    f2 <- paste0("hyper_", name)
    common[[f1]] <- pp$prior
    common[[f2]] <- pp$hyper
  }
  
  # TODO: Change these
  other <- list(
    prior_teff = repvec(c(2, 0, 1), num_uncrt > 0),

    hyper_beta = repvec(c(0.2, 0.2), num_heter > 0),
    hyper_teff = repvec(c(0, 1, 0), num_uncrt > 0),

    teff_obs = repvec(rep(0, num_cases), num_uncrt > 0),
    teff_lb = repvec(rep(0, num_cases), num_uncrt > 0),
    teff_ub = repvec(rep(0, num_cases), num_uncrt > 0)
  )

  # Return
  c(common, other)
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

