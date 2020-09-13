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
  ci <- get_covariate_info(object)
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
  invisible(object)
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
#' @export
#' @inheritParams object_to_model
#' @param x x-axis variable name or a vector to be used as x-axis values
#' @param group_by grouping variable name or a factor to be used for grouping
#' the plot into panels
#' @return a data frame
create_plot_df <- function(object, x = "age", group_by = "id") {
  num_obs <- get_num_obs(object)

  # Get the x-axis variable
  if (is.character(x)) {
    x_name <- x
    x <- get_covariate(object, x)
    check_type(x, "numeric")
  } else {
    x_name <- deparse(substitute(x))
    check_type(x, "numeric")
    check_length(x, num_obs)
  }

  # Get the grouping factor
  if (is.character(group_by)) {
    g_name <- group_by
    group_by <- get_covariate(object, group_by)
    check_type(group_by, "factor")
  } else {
    g_name <- deparse(substitute(group_by))
    check_type(group_by, "factor")
    check_length(group_by, num_obs)
  }

  y <- as.numeric(get_y(object))
  y_name <- get_y_name(object)
  df <- data.frame(group_by, x, y)
  colnames(df) <- c(g_name, x_name, y_name)
  return(df)
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
#'   \item \code{get_covariate_info} returns a list
#'   \item \code{get_component_names} returns an array of names
#'   \item \code{get_num_obs} returns an integer
#'   \item \code{is_f_sampled} returns a boolean value
#'   \item \code{get_obs_model} returns a character string
#'   \item \code{get_ns_covariates} returns a character vector
#'   \item \code{get_y_name} returns the response variable name
#' }
#' @family lgpmodel accessors
NULL

#' @export
#' @param type What types of covariates to include. Must be either
#' "all", "categorical" or "continuous".
#' @rdname model_getters
get_covariate_names <- function(object, type = "all") {
  model <- object_to_model(object)
  allowed <- c("continuous", "categorical", "all")
  idx <- check_allowed(arg = type, allowed = allowed)
  nam1 <- dollar(model@var_names, "x_cont")
  nam2 <- dollar(model@var_names, "x_cat")
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

#' @export
#' @rdname model_getters
get_covariate_info <- function(object) {
  info1 <- get_covariate_info_cont(object)
  info2 <- get_covariate_info_cat(object)
  list(
    continuous = info1,
    categorical = info2
  )
}

#' @export
#' @rdname model_getters
get_component_names <- function(object) {
  comps <- dollar(get_stan_input(object), "components")
  rownames(comps)
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

#' @rdname model_getters
get_y_name <- function(object) {
  model <- object_to_model(object)
  dollar(model@var_names, "y")
}

#' @rdname model_getters
get_covariate_info_cont <- function(object) {
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

#' @rdname model_getters
get_covariate_info_cat <- function(object) {
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
    if (length(a) > 4) {
      a <- "..."
    } else {
      a <- paste(a, collapse = ", ")
    }
    level_names[j] <- a
  }
  df <- data.frame(nam, num_levels, level_names)
  colnames(df) <- c("Factor", "#Levels", "Values")
  rownames(df) <- seq_len(dim(df)[1])
  return(df)
}


#' Functions that access original data stored in a model object
#'
#' @name data_getters
#' @inheritParams object_to_model
#' @param original should the variable be in its original form
#' @param mask_with value to use to mask the originally missing
#' values in the data
#' @return
#' \itemize{
#'   \item \code{get_y} gets response variable measurements
#'   (on their original scaled if \code{original} is \code{TRUE})
#'   \item \code{get_x_cat} gets the categorical covariate matrix
#'   \item \code{get_x_cont} gets the continuous covariate matrix (with
#'   variables on original scale and \code{NaN}s in their original locations
#'   if \code{original} is \code{TRUE})
#'   \item \code{get_covariate} gets a single covariate
#' }
#' @family lgpmodel accessors
NULL

#' @rdname data_getters
get_y <- function(object, original = TRUE) {
  model <- object_to_model(object)
  is_gauss <- get_obs_model(model) == "gaussian"
  nam <- if (is_gauss) "y_cont" else "y_disc"
  y <- dollar(get_stan_input(model), nam)
  scl <- dollar(model@var_scalings, "y")
  rescale <- original && is_gauss
  y <- if (rescale) scl@fun_inv(y) else y
  rownames(y) <- get_y_name(model)
  return(y)
}

#' @rdname data_getters
get_x_cat <- function(object) {
  x <- dollar(get_stan_input(object), "x_cat")
  out <- if (nrow(x) == 0) NULL else x
  return(out)
}

#' @rdname data_getters
get_x_cont <- function(object, original = TRUE, mask_with = NaN) {
  nam <- if (original) "x_cont_unnorm" else "x_cont"
  si <- get_stan_input(object)
  x <- dollar(si, nam)
  if (nrow(x) == 0) {
    out <- NULL
  } else {
    out <- x
    if (original) {
      mask <- as.logical(dollar(si, "x_cont_mask"))
      out[mask] <- mask_with
    }
  }
  return(out)
}

#' @rdname data_getters
#' @param name covariate name
get_covariate <- function(object, name) {
  x_cat <- get_x_cat(object)
  x_cont <- get_x_cont(object, original = TRUE, mask_with = NaN)
  if (name %in% rownames(x_cat)) {
    x <- select_row(x_cat, name)
    return(as.factor(x))
  } else if (name %in% rownames(x_cont)) {
    x <- select_row(x_cont, name)
    return(as.numeric(x))
  } else {
    stop("covariate not found!")
  }
}
