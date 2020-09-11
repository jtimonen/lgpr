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
  str2 <- likelihood_as_str(stan_list$obs_model)
  num_indiv <- length(levels(model@id_variable$values))
  str3 <- paste0(model@id_variable$name, ", ", num_indiv, " levels")
  str4 <- model@time_variable$name
  line1 <- paste0("Formula: ", str1)
  line2 <- paste0("Likelihood: ", str2)
  line3 <- paste0("ID variable: ", str3)
  line4 <- paste0("Time variable: ", str4)
  out <- paste0(line1, "\n", line2, "\n", line3, "\n", line4, "\n")
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
#' @return a data frame
create_plot_df <- function(object) {
  model <- object_to_model(object)
  id_var <- model@id_variable
  time_var <- model@time_variable
  y <- as.numeric(get_y(model))
  y_name <- get_y_name(model)
  df <- data.frame(id_var$values, time_var$values, y)
  colnames(df) <- c(id_var$name, time_var$name, y_name)
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
  comps <- get_stan_input(object)$components
  rownames(comps)
}

#' @export
#' @rdname model_getters
get_num_obs <- function(object) {
  get_stan_input(object)$num_obs
}

#' @export
#' @rdname model_getters
is_f_sampled <- function(object) {
  val <- get_stan_input(object)$is_f_sampled
  as.logical(val)
}

#' @export
#' @rdname model_getters
get_obs_model <- function(object) {
  lh <- get_stan_input(object)$obs_model
  likelihood_as_str(lh)
}

#' @rdname model_getters
get_y_name <- function(object) {
  model <- object_to_model(object)
  model@var_names$y
}

#' @rdname model_getters
get_ns_covariates <- function(object) {
  cinfo <- get_component_info(object)
  x_cont <- get_stan_input(object)$x_cont
  inds <- which(cinfo$ns == 1)
  out <- NULL
  for (idx in inds) {
    icont <- as.numeric(cinfo$cont[idx])
    name <- rownames(x_cont)[icont]
    out <- c(out, name)
  }
  return(out)
}

#' @rdname model_getters
get_covariate_info_cont <- function(object) {
  model <- object_to_model(object)
  vn <- model@var_names
  nam <- vn$x_cont
  if (is.null(nam)) {
    return(NULL)
  }
  num_nan <- rowSums(model@stan_input$x_cont_mask == 1)
  df <- data.frame(nam, num_nan)
  colnames(df) <- c("Variable", "#Missing")
  rownames(df) <- seq_len(dim(df)[1])
  return(df)
}

#' @rdname model_getters
get_covariate_info_cat <- function(object) {
  model <- object_to_model(object)
  vn <- model@var_names
  nam <- vn$x_cat
  if (is.null(nam)) {
    return(NULL)
  }
  num_levels <- model@stan_input$x_cat_num_levels
  levels <- model@var_info$x_cat_levels
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
#' }
#' @family lgpmodel accessors
NULL

#' @rdname data_getters
get_y <- function(object, original = TRUE) {
  model <- object_to_model(object)
  is_gauss <- get_obs_model(model) == "gaussian"
  nam <- if (is_gauss) "y_cont" else "y_disc"
  y <- dollar(get_stan_input(model), nam)
  scl <- model@var_scalings$y
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

#' Get observed disease effect times
#'
#' @export
#' @inheritParams object_to_model
#' @param dra name of the disease-related age variable
#' @return a named vector, where the names are individual IDs
get_teff_obs <- function(object, dra = "diseaseAge") {
  model <- object_to_model(object)
  df_data <- create_plot_df(model)
  ns_names <- get_ns_covariates(model)
  ok <- dra %in% ns_names
  if (!ok) stop(dra, " is not a nonstationary covariate in the model!")
  dis_ages <- select_row(get_x_cont(model), dra)
  ids <- df_data[, 1]
  ages <- df_data[, 2]
  compute_teff_obs(ids, ages, dis_ages)
}
