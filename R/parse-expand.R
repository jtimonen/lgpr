#' Create idx_expand input for Stan
#'
#' @param components the \code{components} input for Stan
#' @param x_cat the \code{x_cat} input for Stan
#' @param x_cont_mask the \code{x_cont_mask} input for Stan
#' @return the \code{idx_expand} input for Stan
create_idx_expand <- function(components, x_cat, x_cont_mask) {
  x_fac <- create_idx_expand_picker(components, x_cat)
  inds <- as.numeric(which(components[, 4] + components[, 7] > 0))
  L <- length(inds)
  if (L == 0) {
    return(x_fac)
  } else {
    i_cont <- components[inds, 9]
    to_red <- x_cont_mask[i_cont, ]
    if (length(i_cont) == 1) {
      to_red <- repvec(to_red, 1)
    }
    idx_mask <- reduce_rows(to_red)
    idx_expand <- create_idx_expand_helper(x_fac, idx_mask)
    return(idx_expand)
  }
}


#' Helper function for creating the idx_expand input for Stan
#'
#' @inheritParams create_idx_expand
#' @return a vector of length \code{n_obs}
create_idx_expand_picker <- function(components, x_cat) {
  n_obs <- dim(x_cat)[2]
  inds <- c(components[, 4], components[, 7])
  inds <- as.numeric(inds[inds != 0])
  J <- length(inds)
  if (J == 0) {
    return(rep(1, n_obs))
  }
  all_same <- all(inds == inds[1])
  if (!all_same) {
    str <- paste(inds, collapse = ", ")
    msg <- paste0(
      "The heter() and uncrt() expressions must have the same ",
      "categorical covariate in every term! ",
      "Found inds = {", str, "}"
    )
    stop(msg)
  }
  x_cat[inds[1], ]
}

#' Check that each row of array is identical
#'
#' @param rows rows
#' @return the indentical row
reduce_rows <- function(rows) {
  nam <- rownames(rows)
  R <- dim(rows)[1]
  r1 <- as.numeric(rows[1, ])
  nam1 <- nam[1]
  n <- length(r1)
  if (R == 1) {
    return(r1)
  }
  for (j in 2:R) {
    rj <- as.numeric(rows[j, ])
    namj <- nam[j]
    s <- sum(rj == r1)
    if (s != n) {
      msg <- paste0(
        "For each term with uncrt() or heter() expressions, ",
        "NaNs of the continuous covariate must be on the same",
        "rows. Found discrepancy between ",
        nam1, " and ", namj, "."
      )
      stop(msg)
    }
  }
  return(r1)
}

#' Helper function for creating the idx_expand input for Stan
#'
#' @inheritParams create_idx_expand
#' @return a vector of length \code{n_obs}
create_idx_expand_helper <- function(x_cat, x_cont_mask) {
  x_cat[x_cont_mask == 1] <- -1
  out <- as.numeric(as.factor(x_cat)) + 1
  out[x_cont_mask == 1] <- 1
  return(out)
}
