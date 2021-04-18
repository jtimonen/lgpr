
#' Graphical posterior or prior predictive checks
#'
#' @export
#' @param data the original data frame
#' @param fun \code{bayesplot} function name
#' @param ... additional arguments passed to the default
#' \code{\link[bayesplot]{pp_check}} method in
#' \code{bayesplot}
#' @param verbose Can this print any messages?
#' @return a \code{ggplot} object
#' @seealso Introduction to graphical posterior predictive checks:
#' \href{here}{https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html}
#'
ppc <- function(fit, data = NULL, fun = default_ppc_fun(fit), verbose = TRUE,
                ...) {
  check_type(fit, "lgpfit")

  # Data is taken from the fit object
  if (!is.null(data)) {
    log_info("The 'data' argument is deprecated and has no effect.")
  }
  data <- get_data(fit)
  check_type(fun, "function")
  y_name <- get_y_name(fit)
  y <- dollar(data, y_name)
  stop("not implemented")
  # y_rep <- get_y_rng(fit, original_scale = TRUE)
  # bayesplot::pp_check(y, y_rep, fun, ...)
}

# Default bayesplot ppc function
default_ppc_fun <- function(object) {
  likelihood <- get_obs_model(object)
  f1 <- bayesplot::ppc_dens_overlay
  f2 <- bayesplot::ppc_hist
  fun <- if (likelihood == "gaussian") f1 else f2
  check_type(fun, "function")
  return(fun)
}
