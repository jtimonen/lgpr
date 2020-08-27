#' Get covariate names
#'
#' @inheritParams object_to_model
#' @param type What types of covariates to include. Must be either
#' "all", "categorical" or "continuous".
#' @return A string where names are separated by a comma.
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

#' Get integer encoding of model components
#'
#' @inheritParams object_to_model
#' @return A data frame.
get_integer_encoding <- function(object) {
  model <- object_to_model(object)
  data.frame(model@stan_input$components)
}

#' Prior summary
#'
#' @inheritParams object_to_model
#' @param digits number of digits to show for floating point numbers
#' @return a data frame
prior_summary <- function(object, digits = 3) {
  model <- object_to_model(object)
  prior_to_df(model@stan_input, digits = digits)
}
