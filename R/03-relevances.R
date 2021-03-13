#' Assess component relevances
#'
#' @description
#' \itemize{
#'   \item \code{relevances} returns a named vector with length equal to
#'   \code{num_comps + 1}
#'   \item \code{relevances.default} is the default method
#'   \item \code{relevances.default_all} returns a \code{data.frame} of size
#'   \code{num_draws} x \code{num_comps + 1}
#' }
#' @param fit an object of class \code{lgpfit}
#' @param reduce a function to apply to reduce the relevances in each
#' draw into one value
#' @param ... currently has no effect
#' @name relevances
NULL

#' @export
#' @rdname relevances
relevances <- function(fit, reduce = function(x) base::mean(x), ...) {
  check_type(fit, "lgpfit")
  relevances.default(fit, reduce, ...)
}

#' @rdname relevances
relevances.default <- function(fit, reduce, ...) {
  df <- relevances.default.all(fit)
  if (is.null(reduce)) {
    return(df)
  } else {
    check_type(reduce, "function")
    df <- apply(df, 2, reduce)
  }
  return(df)
}

#' @rdname relevances
relevances.default.all <- function(fit) {
  pred <- get_pred(fit)
  p_noise <- relevances.default.noise(fit, pred)
  p_comps <- relevances.default.comp(fit, pred)
  num_comps <- ncol(p_comps)
  mult <- t(repvec(1 - p_noise, num_comps))
  p_comps <- mult * p_comps
  df <- cbind(p_comps, p_noise)
  colnames(df) <- c(names(p_comps), "noise")
  return(df)
}

#' Helper functions for relevances.default
#'
#' @description
#' \itemize{
#'   \item \code{relevances.default.noise} computes the noise proportion
#'   in each draws, and returns a vector with length \code{num_draws}
#'   \item \code{relevances.default.comp} divides the signal proportion for
#'   each component in each draw, and returns a data frame with
#'   \code{num_draws} rows and \code{num_comps} colums
#' }
#' @param pred an object of class \linkS4class{FunctionDraws} or
#' \linkS4class{FunctionPosterior}
#' @inheritParams relevances
#' @name relevances_default
#' @return an array or vector
NULL

#' @rdname relevances_default
relevances.default.noise <- function(fit, pred) {
  typ <- class(pred)
  h <- if (typ == "FunctionDraws") pred@h else pred@y_mean
  num_obs <- ncol(h)
  num_draws <- nrow(h)
  y <- get_y(fit, original = TRUE)
  y <- divide_by_num_trials(y, fit)
  resid <- h - repvec(y, num_draws)
  signal_var <- row_vars(h)
  error_var <- rowSums(resid^2) / (num_obs - 1)
  p_signal <- signal_var / (signal_var + error_var)
  return(1 - p_signal)
}

#' @rdname relevances_default
relevances.default.comp <- function(fit, pred) {
  typ <- class(pred)
  f_comp <- if (typ == "FunctionDraws") pred@f_comp else pred@f_comp_mean
  v_list <- lapply(f_comp, row_vars)
  nam <- names(v_list)
  df <- data.frame(v_list)
  colnames(df) <- nam
  df <- normalize_rows(df)
  return(df)
}

#' Select relevant components
#'
#' @description
#' \itemize{
#'   \item \code{select} performs strict selection, returning either \code{TRUE}
#'   or \code{FALSE} for each component.
#'   \item \code{select.integrate} is like \code{select}, but instead of
#'   a fixed threshold, computes probabilistic selection by integrating over
#'   a threshold density.
#'   \item \code{select_freq} performs the selection separately using
#'   each parameter draw and returns the frequency at which each
#'   component was selected.
#'   \item \code{select_freq.integrate} is like \code{select_freq}, but
#'   instead of a fixed threshold, computes probabilistic selection
#'   frequencies by integrating over a threshold density.
#' }
#' @param fit An object of class \code{lgpfit}.
#' @param reduce The \code{reduce} argument for \code{\link{relevances}}.
#' @param threshold Threshold for relevance sum.
#' Must be a value between 0 and 1.
#' @param p A threshold density over interval [0,1].
#' @param h A discretization parameter for computing a quadrature.
#' @param verbose Should this show a progress bar?
#' @param ... Additional arguments to \code{\link{relevances}}.
#' @return See description.
#' @name select
NULL

#' @export
#' @rdname select
select <- function(fit,
                   reduce = function(x) base::mean(x),
                   threshold = 0.95, ...) {
  check_type(fit, "lgpfit")
  check_interval(threshold, 0, 1)
  check_type(reduce, "function")

  r <- relevances(fit, reduce = reduce, ...)
  nam <- names(r)
  L <- length(r)
  r <- array(r, dim = c(1, L))
  r <- as.matrix(r)
  colnames(r) <- nam
  s <- select.all_draws(r, threshold)
  result_df(as.logical(s), "Selected", colnames(s))
}

#' @export
#' @rdname select
select_freq <- function(fit, threshold = 0.95, ...) {
  check_type(fit, "lgpfit")
  check_interval(threshold, 0, 1)
  r <- relevances(fit, reduce = NULL, ...)
  nam <- names(r)
  r <- as.matrix(r)
  colnames(r) <- nam
  sel <- select.all_draws(r, threshold)
  result_df(colMeans(sel), "P(selected)", colnames(sel))
}

#' @export
#' @rdname select
select.integrate <- function(fit,
                             reduce = function(x) base::mean(x),
                             p = function(x) stats::dbeta(x, 100, 5),
                             h = 0.01,
                             verbose = TRUE,
                             ...) {
  selected <- select.all_thresholds(fit, reduce, p, h, verbose, ...)
  prob <- select.compute_integral(selected, p)
  list(
    selected = selected,
    expected = result_df(prob, "P(selected)", colnames(selected))
  )
}

#' @export
#' @rdname select
select_freq.integrate <- function(fit,
                                  p = function(x) stats::dbeta(x, 100, 5),
                                  h = 0.01,
                                  verbose = TRUE,
                                  ...) {
  freqs <- select_freq.all_thresholds(fit, p, h, verbose, ...)

  # Compute integral and return
  prob <- select.compute_integral(freqs, p)
  list(
    freq = freqs,
    expected = result_df(prob, "P(selected)", colnames(freqs))
  )
}


#' Helper functions for component selection
#'
#' @inheritParams select
#' @param rel an array of relevances
#' @name select.helper
NULL


#' @rdname select.helper
select.all_thresholds <- function(fit, reduce, p, h, verbose, ...) {

  # Get relevances
  check_type(reduce, "function")
  check_type(fit, "lgpfit")
  check_type(p, "function")
  check_interval(h, 0, 1)
  check_positive(h)

  rel <- relevances(fit, reduce = reduce, ...)
  D <- length(rel)
  nam <- names(rel)
  rel <- array(rel, dim = c(1, D))
  colnames(rel) <- nam

  # Setup
  H <- seq(0, 1, by = h)
  L <- length(H)
  SEL <- array(0, dim = c(L, D))
  colnames(SEL) <- colnames(rel)
  pb <- progbar_header(L)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  if (verbose) cat(hdr, "\n")

  # Loop
  for (i in 1:L) {
    sel <- select.all_draws(rel, threshold = H[i])
    SEL[i, ] <- sel
    if (verbose) progbar_print(i, idx_print)
  }
  if (verbose) cat("\n")
  return(SEL)
}


#' @rdname select.helper
select_freq.all_thresholds <- function(fit, p, h, verbose, ...) {

  # Get relevances
  check_type(fit, "lgpfit")
  check_type(p, "function")
  check_interval(h, 0, 1)
  check_positive(h)
  rel <- relevances(fit, reduce = NULL, ...)
  D <- ncol(rel)

  # Setup
  H <- seq(0, 1, by = h)
  L <- length(H)
  FREQ <- array(0, dim = c(L, D))
  colnames(FREQ) <- colnames(rel)
  pb <- progbar_header(L)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  if (verbose) cat(hdr, "\n")

  # Loop
  for (i in 1:L) {
    freq <- select.all_draws(rel, threshold = H[i])
    FREQ[i, ] <- colMeans(freq)
    if (verbose) progbar_print(i, idx_print)
  }
  if (verbose) cat("\n")
  return(FREQ)
}


#' @rdname select.helper
#' @param arr an array of frequencies or selections
select.compute_integral <- function(arr, p) {
  L <- nrow(arr)
  D <- ncol(arr)
  H <- seq(0, 1, length.out = L)
  P <- p(H)
  P <- P / sum(P)
  E <- t(repvec(P, D))
  colSums(E * arr)
}

#' @rdname select.helper
select.all_draws <- function(rel, threshold) {
  S <- nrow(rel)
  J <- ncol(rel)
  sel <- array(0, dim = c(S, J))
  for (s in seq_len(S)) {
    inds <- select.one_draw(rel[s, ], threshold)
    sel[s, inds] <- 1
  }
  colnames(sel) <- colnames(rel)
  return(sel)
}

#' @rdname select.helper
select.one_draw <- function(rel, threshold) {
  J <- length(rel)
  p_noise <- rel[J]
  rel <- as.numeric(rel[1:(J - 1)])
  s <- sort(rel, decreasing = TRUE, index.return = TRUE)
  rel <- s$x
  s_ix <- s$ix
  if (p_noise >= threshold) {
    return(J)
  }
  for (j in seq_len(J - 1)) {
    h <- p_noise + sum(rel[1:j])
    if (h >= threshold) {
      i_sel <- c(s_ix[1:j], J)
      return(i_sel)
    }
  }
  return(1:J)
}

#' @rdname select.helper
#' @param val values
#' @param val_name name of variable that \code{val} presents
#' @param comp_names component names
result_df <- function(val, val_name, comp_names) {
  df <- data.frame(val)
  colnames(df) <- c("Component")
  rownames(df) <- comp_names
  return(df)
}
