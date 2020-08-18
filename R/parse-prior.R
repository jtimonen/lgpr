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
    prior_alpha = repvec(c(3,0), num_comps),
    prior_ell = repvec(c(6,0), num_ell),
    prior_wrp = repvec(c(6,0), num_ns),
    prior_sigma = repvec(c(6,0), num_sigma),
    prior_phi = repvec(c(6,0), num_phi),
    prior_teff = repvec(c(2,0,1), num_uncrt > 0),
    
    hyper_alpha = repvec(c(20,1,0), num_comps),
    hyper_ell = repvec(c(1,1,0), num_ell),
    hyper_wrp = repvec(c(1,1,0), num_ns),
    hyper_sigma = repvec(c(1,1,0), num_sigma),
    hyper_phi = repvec(c(1,1,0), num_phi),
    hyper_beta = repvec(c(0.2, 0.2), num_heter > 0),
    hyper_teff = repvec(c(0,1,0), num_uncrt > 0),
    
    teff_obs = repvec(rep(0, num_cases), num_uncrt > 0),
    teff_lb = repvec(rep(0, num_cases), num_uncrt > 0),
    teff_ub = repvec(rep(0, num_cases), num_uncrt > 0)
  )
}

