#' Print a model summary
#'
#' @export
#' @inheritParams object_to_model
#' @return
#' \itemize{
#'   \item \code{\link{model_summary_brief}} returns a string
#'   \item \code{\link{model_summary}} prints the summary and
#' returns \code{object} invisibly
#'   \item \code{\link{print_stan_input}} prints the stan input
#'   and returns \code{object} invisibly
#' }
model_summary <- function(object) {
  model <- object_to_model(object)
  brief <- model_summary_brief(object)
  cat(brief)
  cat("\n")
  print(get_component_info(object))
  cat("\n")
  print(param_summary(model))
  invisible(object)
}

#' @export
#' @rdname model_summary
model_summary_brief <- function(object) {
  model <- object_to_model(object)
  stan_list <- get_stan_input(object)
  str1 <- as.character(model@model_formula)
  str2 <- likelihood_as_str(stan_list$obs_model)
  line1 <- paste0("Formula: ", str1)
  line2 <- paste0("Likelihood: ", str2)
  out <- paste0(line1, "\n", line2, "\n")
  return(out)
}

#' @export
#' @rdname model_summary
print_stan_input <- function(object) {
  alist <- get_stan_input(object)
  print_list(alist)
  invisible(model)
}

#' Parameter summary (priors etc.)
#'
#' @export
#' @inheritParams object_to_model
#' @param digits number of digits to show for floating point numbers
#' @return data frame
param_summary <- function(object, digits = 3) {
  model <- object_to_model(object)
  prior_to_df(model@stan_input, digits = digits)
}

#' @export
#' @rdname param_summary
prior_summary <- function(object, digits = 3) {
  param_summary(object, digits)
}

#' Functions that access model properties
#'
#' @name model_getters
#' @inheritParams object_to_model
#' @return
#' \itemize{
#'   \item \code{get_covariate_names} returns a string where names
#'   are separated by a comma
#'   \item \code{get_stan_input} returns a list
#'   \item \code{get_component_info} returns a data frame
#'   \item \code{get_component_names} returns an array of names
#'   \item \code{get_num_obs} returns an integer
#' }
NULL

#' @export
#' @param type What types of covariates to include. Must be either
#' "all", "categorical" or "continuous".
#' @rdname model_getters
get_covariate_names <- function(object, type = "all") {
  model <- object_to_model(object)
  allowed <- c("continuous", "categorical", "all")
  idx <- check_allowed(arg = type, allowed = allowed)
  nam1 <- model@var_names$x_cont
  nam2 <- model@var_names$x_cat
  nam3 <- c(nam1, nam2)
  names <- list(nam1, nam2, nam3)
  out <- names[[idx]]
  paste(out, collapse = ", ")
}

#' @export
#' @rdname model_getters
get_stan_input <- function(object) {
  model <- object_to_model(object)
  return(model@stan_input)
}

#' @export
#' @rdname model_getters
get_component_info <- function(object) {
  comps <- get_stan_input(object)$components
  a <- cbind(comps[, 1:2], comps[, 4:9])
  Component <- rownames(a)
  a <- cbind(Component, a)
  rownames(a) <- NULL
  data.frame(a)
}

#' @export
#' @rdname model_getters
get_component_names <- function(object) {
  comps <- get_stan_input(object)$components
  rownames(comps)
}

#' @export
#' @rdname model_getters
get_num_obs <- function(object) {
  get_stan_input(object)$num_obs
}
