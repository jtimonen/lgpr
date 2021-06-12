#' @export
#' @describeIn lgpmodel Print information and summary about the object.
#' Returns \code{object} invisibly.
setMethod("show", "lgpmodel", function(object) {
  desc <- class_info("lgpmodel")
  cat(desc)
  cat("\n")
  model_summary(object)
})

#' @export
#' @describeIn lgpmodel Get a parameter summary (bounds and
#' priors). Returns a \code{data.frame}.
#' @param digits number of digits to show for floating point numbers
setMethod("parameter_info", "lgpmodel", function(object, digits = 3) {
  si <- get_stan_input(object)
  prior_to_df(si, digits = digits)
})

#' @export
#' @describeIn lgpmodel Get a data frame with information about each model
#' component.
setMethod("component_info", "lgpmodel", function(object) {
  comps <- get_component_encoding(object)
  nams <- colnames(comps)
  p1 <- ensure_2dim(comps[, 1:2])
  p2 <- ensure_2dim(comps[, 4:9])
  a <- cbind(p1, p2)
  Component <- rownames(comps)
  a <- cbind(Component, a)
  colnames(a) <- c("Component", nams[1:2], nams[4:9])
  rownames(a) <- NULL
  data.frame(a)
})

#' @export
#' @describeIn lgpmodel Get covariate information.
setMethod("covariate_info", "lgpmodel", function(object) {
  info1 <- covariate_info.cont(object)
  info2 <- covariate_info.cat(object)
  list(continuous = info1, categorical = info2)
})

#' @export
#' @describeIn lgpmodel Get names of model components.
setMethod("component_names", "lgpmodel", function(object) {
  rownames(get_component_encoding(object))
})

#' @export
#' @describeIn lgpmodel Determine if inference of the model requires sampling
#' the latent signal \code{f} (and its components).
setMethod("is_f_sampled", "lgpmodel", function(object) {
  object@sample_f
})

#' Print a model summary.
#'
#' @export
#' @param object a model or fit
#' @param digits number of digits to round floats to
#' @return \code{object} invisibly.
model_summary <- function(object, digits = 3) {
  model <- object_to_model(object)

  # Helper function
  model_info <- function(object) {
    model <- object_to_model(object)
    stan_list <- get_stan_input(object)
    str1 <- as.character(model@model_formula)
    str2 <- likelihood_as_str(dollar(stan_list, "obs_model"))
    dat <- get_data(object)
    N <- nrow(dat)
    D <- ncol(dat)
    line1 <- paste0("Formula: ", str1)
    line2 <- paste0("Likelihood: ", str2)
    line3 <- paste0("Data: ", N, " observations, ", D, " variables")
    out <- paste0(line1, "\n", line2, "\n", line3, "\n")
    return(out)
  }

  brief <- model_info(model)
  cat(brief)
  cat("\n")
  print(component_info(model))
  cat("\n")
  ci <- covariate_info(model)
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
  print(parameter_info(model, digits))
  bi <- beta_teff_idx_info(model)
  if (!is.null(bi)) {
    cat("\n")
    print(bi)
  }
  cat("\n")
  cat(misc_info(model), "\n")
  invisible(object)
}

#' @export
#' @rdname model_summary
param_summary <- function(object, digits = 3) {
  model <- object_to_model(object)
  parameter_info(model, digits)
}

beta_teff_idx_info <- function(object) {
  model <- object_to_model(object)
  a <- dollar(model@info, "caseid_map")
  if (is.null(a)) {
    return(a)
  }
  a <- t(a)
  rownames(a)[2] <- "beta_or_teff_param_idx"
  a <- data.frame(a)
  colnames(a) <- NULL
  return(a)
}

misc_info <- function(object) {
  model <- object_to_model(object)
  info <- model@info
  desc <- paste0(
    "Created on ", dollar(info, "created"), " with lgpr ",
    dollar(info, "lgpr_version"), "."
  )
  return(desc)
}

# Categorial covariate information
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

# Continuous covariate information
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

# Get the c_chat Stan input or a vector of zeros
get_chat <- function(object) {
  si <- get_stan_input(object)
  model <- object_to_model(object)
  if (is_f_sampled(model)) {
    c_hat <- dollar(si, "c_hat")
  } else {
    N <- get_num_obs(model)
    c_hat <- rep(0.0, N)
  }
  return(c_hat)
}

# Get response variable measurements on original or normalized scale
get_y <- function(object, original = TRUE) {
  if (original) {
    y_name <- get_y_name(object)
    dat <- get_data(object)
    return(dollar(dat, y_name))
  }
  if (is_f_sampled(object)) {
    stop(
      "Response variable is not normalized if f is sampled! Set ",
      "original = TRUE."
    )
  }
  si <- get_stan_input(object)
  out <- as.vector(dollar(si, "y_norm"))
  return(out)
}

# Get response variable name
get_y_name <- function(object) {
  model <- object_to_model(object)
  dollar(model@var_names, "y")
}

# Get the Stan model used by a model
get_stan_model <- function(object) {
  model <- object_to_model(object)
  if (is_f_sampled(model)) {
    model_name <- "lgp_latent"
  } else {
    model_name <- "lgp"
  }
  stanmodels[[model_name]] # global variable (list of all pkg models)
}

# Get Stan input
get_stan_input <- function(object) {
  model <- object_to_model(object)
  return(model@stan_input)
}

# Get integer matrix encoding component types
get_component_encoding <- function(object) {
  si <- get_stan_input(object)
  dollar(si, "components")
}

# Get raw original data
get_data <- function(object) {
  model <- object_to_model(object)
  return(model@data)
}

# Get number of observations
get_num_obs <- function(object) {
  dollar(get_stan_input(object), "num_obs")
}

# Get number of components
get_num_comps <- function(object) {
  dollar(get_stan_input(object), "num_comps")
}

# Get observation model (human readable string)
get_obs_model <- function(object) {
  lh <- dollar(get_stan_input(object), "obs_model")
  likelihood_as_str(lh)
}

# Get number of trials (binomial or BB model)
get_num_trials <- function(object) {
  num_trials <- dollar(get_stan_input(object), "y_num_trials")
  as.vector(num_trials)
}


#' Helper function for plots
#'
#' @param object model or fit
#' @param x x-axis variable name
#' @param group_by grouping variable name (use \code{NULL} for no grouping)
#' @return a data frame
create_plot_df <- function(object, x = "age", group_by = "id") {

  # Get x-axis variable
  dat <- get_data(object)
  x_name <- x
  x <- dollar(dat, x_name)
  check_type(x, "numeric")

  # Get grouping factor
  x_grp <- create_grouping_factor(dat, group_by) # util

  # Get response
  y <- get_y(object, original = TRUE)
  y_name <- get_y_name(object)
  df <- data.frame(x_grp, x, y)
  group_by <- if (is.na(group_by)) "group__" else group_by
  colnames(df) <- c(group_by, x_name, y_name)
  return(df)
}
