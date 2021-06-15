#' Graphical posterior predictive checks
#'
#' @export
#' @param fit An object of class \linkS4class{lgpfit} that can been created
#' with \code{sample_f=TRUE}.
#' @param data the original data frame (deprecated argument with no
#' effect, now obtained from fit object)
#' @param fun \code{bayesplot} function name
#' @param ... additional arguments passed to the default
#' \code{\link[bayesplot]{pp_check}} method in
#' \code{bayesplot}
#' @param verbose Can this print any messages?
#' @return a \code{ggplot} object
#' @seealso Introduction to graphical posterior predictive checks:
#' \href{https://CRAN.R-project.org/package=bayesplot/vignettes/graphical-ppcs.html}{here}.
#' Prior predictive check can be done by calling
#' \code{\link{prior_pred}} and then \code{bayesplot::pp_check()}.
ppc <- function(fit, data = NULL, fun = default_ppc_fun(fit), verbose = TRUE,
                ...) {
  check_type(fit, "lgpfit")
  if (!is_f_sampled(fit)) {
    stop("fit has been created with sample_f = FALSE")
  }

  # Data is taken from the fit object
  if (!is.null(data)) {
    log_info("The 'data' argument is deprecated and has no effect.")
  }
  data <- get_data(fit)
  check_type(fun, "function")
  y_name <- get_y_name(fit)
  y <- dollar(data, y_name)
  y_rep <- draw_pred(fit)
  bayesplot::pp_check(y, y_rep, fun, ...)
}

# Default bayesplot ppc function
default_ppc_fun <- function(object) {
  bayesplot::ppc_dens_overlay
}
