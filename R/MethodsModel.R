#' Print a model summary
#'
#' @export
#' @inheritParams object_to_model
#' @return
#' \itemize{
#'   \item \code{\link{model_summary_brief}} returns a string
#'   \item \code{\link{model_summary}} prints the summary and
#' returns \code{object} invisibly
#'   and returns \code{object} invisibly
#' }
model_summary <- function(object) {
  model <- object_to_model(object)
  brief <- model_summary_brief(object)
  cat(brief)
  cat("\n")
  print(component_info(object))
  cat("\n")
  ci <- covariate_info(object)
  info_cont <- dollar(ci, "continuous")
  info_cat <- dollar(ci, "categorical")
  if (!is.null(info_cont)) {
    print(info_cont)
    cat("\n")
  }
  if (!is.null(info_cat)) {
    print(info_cat)
    cat("\n")
  }
  print(param_summary(model))
  invisible(object)
}

#' @export
#' @rdname model_summary
model_summary_brief <- function(object) {
  model <- object_to_model(object)
  stan_list <- get_stan_input(object)
  str1 <- as.character(model@model_formula)
  str2 <- likelihood_as_str(dollar(stan_list, "obs_model"))
  dat <- get_data(object)
  N <- nrow(dat)
  D <- ncol(dat)
  line1 <- paste0("Formula: ", str1)
  line2 <- paste0("Likelihood: ", str2)
  line3 <- paste("Data:", N, "observations,", D, " variables")
  out <- paste0(line1, "\n", line2, "\n", line3, "\n")
  return(out)
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

#' Helper function for plots
#'
#' @inheritParams object_to_model
#' @param x x-axis variable name
#' @param group_by grouping variable name
#' @return a data frame
create_plot_df <- function(object, x = "age", group_by = "id") {

  # Get x-axis variable
  dat <- get_data(object)
  x_name <- x
  x <- dollar(dat, x_name)
  check_type(x, "numeric")

  # Get grouping factor
  g_name <- group_by
  group_by <- dollar(dat, g_name)
  check_type(group_by, "factor")

  # Get response
  y <- get_y(object, original = TRUE)
  y_name <- get_y.name(object)
  df <- data.frame(group_by, x, y)
  colnames(df) <- c(g_name, x_name, y_name)
  return(df)
}

#' Information about covariates used in a model
#'
#' @inheritParams object_to_model
#' @return a list or a string
#' @name covariate_info
NULL

#' @rdname covariate_info
covariate_info <- function(object) {
  info1 <- covariate_info.cont(object)
  info2 <- covariate_info.cat(object)
  list(
    continuous = info1,
    categorical = info2
  )
}

#' @rdname covariate_info
covariate_info.cat <- function(object) {
  model <- object_to_model(object)
  vn <- model@var_names
  nam <- dollar(vn, "x_cat")
  if (is.null(nam)) {
    return(NULL)
  }
  num_levels <- dollar(model@stan_input, "x_cat_num_levels")
  levels <- dollar(model@var_info, "x_cat_levels")
  level_names <- c()
  J <- length(nam)
  for (j in seq_len(J)) {
    a <- levels[[j]]
    a <- if (length(a) > 4) "..." else paste(a, collapse = ", ")
    level_names[j] <- a
  }
  df <- data.frame(nam, num_levels, level_names)
  colnames(df) <- c("Factor", "#Levels", "Values")
  rownames(df) <- seq_len(dim(df)[1])
  return(df)
}

#' @rdname covariate_info
covariate_info.cont <- function(object) {
  model <- object_to_model(object)
  vn <- model@var_names
  nam <- dollar(vn, "x_cont")
  if (is.null(nam)) {
    return(NULL)
  }
  mask <- dollar(model@stan_input, "x_cont_mask")
  num_nan <- rowSums(mask == 1)
  df <- data.frame(nam, num_nan)
  colnames(df) <- c("Variable", "#Missing")
  rownames(df) <- seq_len(dim(df)[1])
  return(df)
}


#' Information about the model components
#'
#' @inheritParams object_to_model
#' @return a list or a string
#' @name component_info
NULL

#' @export
#' @rdname component_info
component_info <- function(object) {
  comps <- dollar(get_stan_input(object), "components")
  nams <- colnames(comps)
  p1 <- ensure_2dim(comps[, 1:2])
  p2 <- ensure_2dim(comps[, 4:9])
  a <- cbind(p1, p2)
  Component <- rownames(comps)
  a <- cbind(Component, a)
  colnames(a) <- c("Component", nams[1:2], nams[4:9])
  rownames(a) <- NULL
  data.frame(a)
}

#' @rdname component_info
component_names <- function(object) {
  comps <- dollar(get_stan_input(object), "components")
  rownames(comps)
}


#' Get response variable measurements
#'
#' @inheritParams object_to_model
#' @param original should the measuments be on their original scale?
#' @name get_y
#' @return a vector
NULL

#' @rdname get_y
get_y <- function(object, original = TRUE) {
  if (original) {
    y_name <- get_y.name(object)
    dat <- get_data(object)
    return(dollar(dat, y_name))
  }
  model <- object_to_model(object)
  is_gauss <- get_obs_model(model) == "gaussian"
  nam <- if (is_gauss) "y_cont" else "y_disc"
  si <- get_stan_input(object)
  out <- as.vector(dollar(si, nam))
  return(out)
}

#' @rdname get_y
get_y.name <- function(object) {
  model <- object_to_model(object)
  dollar(model@var_names, "y")
}


#' Functions that access model properties
#'
#' @name model_getters
#' @inheritParams object_to_model
#' @return
#' \itemize{
#'   \item \code{get_stan_input} returns a list
#'   \item \code{get_num_obs} returns the number of observations
#'   \item \code{get_obs_model} returns the obs. model as a string
#'   \item \code{get_y_name} returns the response variable name
#'   \item \code{get_data} returns the original unmodified data frame
#'   \item \code{is_f_sampled} returns a boolean value
#' }
NULL

#' @export
#' @rdname model_getters
get_stan_input <- function(object) {
  model <- object_to_model(object)
  return(model@stan_input)
}

#' @export
#' @rdname model_getters
get_data <- function(object) {
  model <- object_to_model(object)
  return(model@data)
}

#' @export
#' @rdname model_getters
get_num_obs <- function(object) {
  dollar(get_stan_input(object), "num_obs")
}

#' @export
#' @rdname model_getters
is_f_sampled <- function(object) {
  val <- dollar(get_stan_input(object), "is_f_sampled")
  as.logical(val)
}

#' @export
#' @rdname model_getters
get_obs_model <- function(object) {
  lh <- dollar(get_stan_input(object), "obs_model")
  likelihood_as_str(lh)
}
