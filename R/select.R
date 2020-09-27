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

