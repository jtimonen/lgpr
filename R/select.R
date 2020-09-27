#' Select relevant components
#'
#' @export
#' @param fit An object of class \code{lgpfit}.
#' @param reduce The \code{reduce} argument for \code{\link{relevances}}.
#' @param threshold Threshold for relevance sum.
#' Must be a value between 0 and 1.
#' @param ... Additional arguments to \code{\link{relevances}}.
#' @return A binary array of shape \code{S} x \code{num_comps+1},
#' where 1 denotes selected components and 0 not selected. \code{S} is
#' equal to \code{num_draws} if \code{reduce} is NULL and otherwise one.
select <- function(fit, reduce = mean, threshold = 0.95, ...) {
  check_type(fit, "lgpfit")
  check_interval(threshold, 0, 1)

  r <- relevances(fit, reduce = reduce, ...)
  nam <- names(r)
  if (!is.null(reduce)) {
    L <- length(r)
    r <- array(r, dim = c(1, L))
  }
  r <- as.matrix(r)
  colnames(r) <- nam
  sel <- select.all_draws(r, threshold)
  return(sel)
}

#' Selection probabilities using a fixed threshold
#'
#' @param rel An array of shape \code{num_draws} x \code{num_comps+1},
#' where rows sum to 1.
#' @inheritParams select
#' @return a binary array with same shape as \code{rel}
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

#' Select relevant components
#'
#' @param rel A vector of component relevances that sums to 1, noise being
#' the last element.
#' @inheritParams select
#' @return indices of selected components (including "noise" always)
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

#' Compute expected selection frequencies by integrating the threshold
#' over [0,1]
#' @export
#' @param object An object of class \code{lgpfit}.
#' @param p The threshold density over interval [0,1].
#' @param h A discretization parameter for computing a quadrature.
#' @param show_progbar Should this show a progress bar?
#' @param ... additional arguments to \code{\link{relevances}}
#' @return A named list
#' @name select_integrate
NULL

#' @rdname select_integrate
select_integrate <- function(fit,
                        p = function(x) {
                          stats::dbeta(x, 100, 5)
                        },
                        h = 0.01,
                        show_progbar = TRUE,
                        ...) {
  check_type(fit, "lgpfit")
  check_type(p, "function")
  check_interval(h, 0, 1)
  check_positive(h)

  # Compute selection frequencies for all threshold values
  freqs <- select_integrate.all(fit, p, h, show_progbar, ...)

  # Compute integral and return
  ef <- select_integrate.integrate(freqs, p)
  list(
    freqs = freqs,
    expected_freqs = dollar(ef, "expectations"),
    threshold_dens = dollar(ef, "threshold_dens")
  )
}

#' @rdname select_integrate
select_integrate.all <- function(fit, p, h, show_progbar, ...) {
  check_type(fit, "lgpfit")
  check_type(p, "function")
  check_interval(h, 0, 1)
  check_positive(h)

  # Get relevances
  rel <- relevances(fit, reduce = NULL, ...)
  D <- ncol(rel)

  # Setup
  H <- seq(0, 1, by = h)
  L <- length(H)
  freqs <- array(0, dim = c(L, D))
  colnames(freqs) <- colnames(rel)
  pb <- progbar_header(L, width = 4)
  hdr <- dollar(pb, "header")
  idx_print <- dollar(pb, "idx_print")
  if (show_progbar) cat(hdr, "\n")

  # Loop
  for (i in 1:L) {
    sel <- select(fit, reduce = NULL, threshold = H[i])
    freqs[i, ] <- colMeans(sel)
    if (show_progbar) {
      if (i %in% idx_print) cat("=")
      if (i == L) cat("\n")
    }
  }
  return(freqs)
}

#' @rdname select_integrate
#' @param freqs an object returned by \code{select_integrate.all}
select_integrate.integrate <- function(freqs, p) {
  L <- nrow(freqs)
  D <- ncol(freqs)
  H <- seq(0, 1, length.out = L)
  P <- p(H)
  list(
    threshold_dens = P,
    expectations = H
  )
}
