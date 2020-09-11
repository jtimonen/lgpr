#' Assess component relevances
#'
#' @description
#' \itemize{
#'   \item \code{relevances} assesses covariate relevances
#'   \item \code{relevances.default} is the only currently implemented method
#' }
#' @param fit an object of class \code{lgpfit}
#' @name relevances
NULL

#' @export
#' @rdname relevances
relevances <- function(fit) {
  check_type(fit, "lgpfit")
  relevances.default(fit)
}

#' @rdname relevances
relevances.default <- function(fit) {
  p_noise <- p_explained_noise(fit)
  p_noise <- mean(p_noise)
  p_comps <- p_explained_components(fit)
  p_comps <- (1 - p_noise) * p_comps
  Component <- c(names(p_comps), "noise")
  Relevance <- c(as.numeric(p_comps), p_noise)
  data.frame(Component, Relevance)
}

#' Determining the variance explained by noise and the signal components
#'
#' @description
#' \itemize{
#'   \item \code{residuals} computes the residuals in each sample
#'   \item \code{p_explained_noise} computes the noise proportion in each sample
#'   \item \code{p_explained_components} divides the signal proportion for
#'   each component
#'   \item \code{component_variances} computes variances in each
#'   sample for each component and returns a \code{data.frame}
#' }
#' @inheritParams relevances
#' @name p_explained
#' @return an array or vector
NULL

#' @rdname p_explained
residuals <- function(fit) {
  y <- get_y(fit, original = FALSE)
  g_total <- get_g(fit)
  num_draws <- nrow(g_total)
  g_total - repvec(y, num_draws)
}

#' @rdname p_explained
p_explained_noise <- function(fit) {
  g_total <- get_g(fit)
  num_obs <- ncol(g_total)
  resid <- residuals(fit)
  signal_var <- row_vars(g_total)
  error_var <- rowSums(resid^2) / (num_obs - 1)
  p_signal <- signal_var / (signal_var + error_var)
  return(1 - p_signal)
}

#' @rdname p_explained
p_explained_components <- function(fit) {
  df <- component_variances(fit)
  r <- colMeans(df)
  r / sum(r)
}

#' @rdname p_explained
component_variances <- function(fit) {
  flag <- is_f_sampled(fit)
  field <- if (flag) NULL else "mean"
  f_comps <- get_f_components(fit, field)
  v_list <- lapply(f_comps, row_vars)
  nam <- names(v_list)
  df <- data.frame(v_list)
  colnames(df) <- nam
  return(df)
}
