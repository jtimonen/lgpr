#' Assess component relevances
#'
#' \itemize{
#'   \item \code{relevances} assesses covariate relevances
#'   \item \code{relevances.f_marginalized} assesses covariate relevances for a
#'   model where the latent function \code{f} was not sampled
#'   \item \code{relevances.f_sampled} assesses covariate relevances for a
#'   model where the latent function \code{f} was sampled
#' }
#' @param fit an object of class \code{lgpfit}
#' @name relevances
NULL

#' @rdname relevances
relevances <- function(fit) {
  check_type(fit, "lgpfit")
  flag <- is_f_sampled(fit)
  if (flag) {
    out <- relevances.f_sampled(fit)
  } else {
    out <- relevances.f_marginalized(fit)
  }
  return(out)
}

#' @rdname relevances
relevances.f_sampled <- function(fit) {
  get_component_info(fit)
}

#' @rdname relevances
relevances.f_marginalized <- function(fit) {
  f <- get_f(fit)
  f_post <- dollar(f, "f")
  f_post
}
