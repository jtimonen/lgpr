#' Print a model summary
#'
#' @export
#' @inheritParams object_to_model
#' @return
#' \itemize{
#'   \item \code{\link{model_summary.brief}} returns a string
#'   \item \code{\link{model_summary}} prints the summary and
#' returns \code{object} invisibly
#'   and returns \code{object} invisibly
#' }
model_summary <- function(object) {
  model <- object_to_model(object)
  brief <- model_summary.brief(object)
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

#' @rdname model_summary
model_summary.brief <- function(object) {
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
  y_name <- get_y_name(object)
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

#' Get the used GP mean vector
#'
#' @inheritParams object_to_model
#' @return a vector with length equal to number of observations
get_chat <- function(object) {
  si <- get_stan_input(object)
  dollar(si, "c_hat")
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
    y_name <- get_y_name(object)
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
get_y_name <- function(object) {
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
#'   \item \code{get_num_trials} returns the vector of numbers of trials
#'   \item \code{is_f_sampled} returns a boolean value
#' }
NULL

#' @rdname model_getters
get_stan_input <- function(object) {
  model <- object_to_model(object)
  return(model@stan_input)
}

#' @rdname model_getters
get_data <- function(object) {
  model <- object_to_model(object)
  return(model@data)
}

#' @rdname model_getters
get_num_obs <- function(object) {
  dollar(get_stan_input(object), "num_obs")
}

#' @rdname model_getters
is_f_sampled <- function(object) {
  val <- dollar(get_stan_input(object), "is_f_sampled")
  as.logical(val)
}

#' @rdname model_getters
get_obs_model <- function(object) {
  lh <- dollar(get_stan_input(object), "obs_model")
  likelihood_as_str(lh)
}

#' @rdname model_getters
get_num_trials <- function(object) {
  num_trials <- dollar(get_stan_input(object), "y_num_trials")
  as.vector(num_trials)
}

#' Compute constant kernel matrices which don't depend on parameters
#'
#' @description
#' \itemize{
#'    \item \code{const_kernels} computes all constant kernels
#'    \item \code{const_kernels.decompositions} computes their exact
#'    rank decompositions
#' }
#' @inheritParams object_to_model
#' @param STREAM external pointer
#' @name const_kernels
#' @return a list of matrices
NULL

#' @rdname const_kernels
const_kernels <- function(object, STREAM = get_stream()) {
  si <- get_stan_input(object)
  n <- get_num_obs(object)
  x_cat <- matrix_to_list(dollar(si, "x_cat"))
  x_cont_mask <- matrix_to_list(dollar(si, "x_cont_mask"))
  x_cat_num_levels <- dollar(si, "x_cat_num_levels")
  components <- matrix_to_list(dollar(si, "components"))
  K_const <- STAN_kernel_const_all(
    n, n,
    x_cat, x_cat,
    x_cont_mask, x_cont_mask,
    x_cat_num_levels, components,
    STREAM
  )
  return(K_const)
}

#' @rdname const_kernels
const_kernels.decompose <- function(object, STREAM = get_stream()) {
  K_const <- const_kernels(object, STREAM)
  si <- get_stan_input(object)
  x_cat_num_levels <- dollar(si, "x_cat_num_levels")
  components <- matrix_to_list(dollar(si, "components"))
  ranks <- STAN_ranks(components, x_cat_num_levels, STREAM)
  list(
    ranks = ranks,
    Delta = STAN_delta_matrix(K_const, ranks, STREAM),
    Theta = STAN_theta_matrix(K_const, ranks, STREAM)
  )
}

#' @rdname const_kernels
const_kernels.decompositions <- function(object, STREAM = get_stream()) {
  rank_dec <- const_kernels.decompose(object, STREAM)
  ranks <- dollar(rank_dec, "ranks")
  Delta <- dollar(rank_dec, "Delta")
  Theta <- dollar(rank_dec, "Theta")
  num_comps <- length(ranks)
  decompositions <- list()
  idx <- 1
  for (j in seq_len(num_comps)) {
    r <- ranks[j]
    if (r > 0) {
      delta_diag <- Delta[idx:(idx + r - 1)]
      if (length(delta_diag) == 1) delta_diag <- as.matrix(delta_diag)
      Delta_j <- diag(delta_diag)
      Theta_j <- Theta[, idx:(idx + r - 1), drop = FALSE]
      decompositions[[j]] <- list(Delta = Delta_j, Theta = Theta_j)
      idx <- idx + r
    } else {
      decompositions[[j]] <- list(Delta = NULL, Theta = NULL)
    }
  }
  return(decompositions)
}


#' @rdname const_kernels
#' @param decompositions a list returned by \code{const_kernel.decompositions}
const_kernels.reconstruct <- function(decompositions, STREAM = get_stream()) {
  num_comps <- length(decompositions)
  K_rec <- list()
  for (j in seq_len(num_comps)) {
    D <- dollar(decompositions[[j]], "Delta")
    TH <- dollar(decompositions[[j]], "Theta")
    if (is.null(D)) {
      K_rec[[j]] <- NA
    } else {
      K_rec[[j]] <- TH %*% D %*% t(TH)
    }
  }
  return(K_rec)
}
