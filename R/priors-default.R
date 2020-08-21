#' Default prior
#'
#' @param name parameter name
#' @return a named list
prior_default <- function(name) {
  
  if (name == 'alpha') {
    prior <- student_t(20)
  } else if (name %in% c('ell', 'sigma', 'phi')) {
    prior <- log_normal(0, 1)
  } else if (name == 'wrp') {
    prior <- log_normal(1, 1)
  } else if (name == 'beta') {
    prior <- beta(a = 0.2, b = 0.2)
  } else if (name == 'effect_time') {
    prior <- 'effect_time_prior'
  } else {
    stop("invalid parameter name '", name, "'")
  }
  return(prior)

}