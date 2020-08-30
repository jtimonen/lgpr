#' Parse the given modeling options
#'
#' @param options A named list with the following possible fields:
#' \itemize{
#'   \item \code{sample_f} Determines if the function values are be sampled
#'   (must be \code{TRUE} if likelihood is not \code{"gaussian"}).
#'   \item \code{skip_generated} If this is true, the generated quantities
#'   block of Stan is skipped.
#'   \item \code{delta} Amount of added jitter to ensure positive definite
#'   covariance matrices.
#' }
#' @return a named list of parsed options
parse_options <- function(options = NULL) {
  input <- options

  # Set defaults
  opts <- list(
    skip_generated = FALSE,
    sample_f = FALSE,
    delta = 1e-8
  )

  # Replace defaults if found from input
  for (opt_name in names(opts)) {
    if (opt_name %in% names(input)) {
      opts[[opt_name]] <- input[[opt_name]]
    }
  }

  # Format for Stan input
  list(
    is_generated_skipped = as.numeric(opts$skip_generated),
    is_f_sampled = as.numeric(opts$sample_f),
    delta = opts$delta
  )
}


#' Check that data contains a variable with a certain name
#'
#' @param var_name the variable to be searched for
#' @param data an object of class \code{data.frame}
#' @return \code{TRUE} if the variable is found
check_in_data <- function(var_name, data) {
  d_names <- colnames(data)
  ok <- (var_name %in% d_names)
  if (!ok) {
    str <- paste(d_names, collapse = ", ")
    msg <- paste0(
      "The variable '", var_name, "' not found in <data>! ",
      " Found data columns = {", str, "}."
    )
    stop(msg)
  }
  return(TRUE)
}

#' Create the function that does a standardizing transform and its inverse
#'
#' @param y variable measurements
#' @param name variable name
#' @return an object of class \linkS4class{lgpscaling}
create_scaling <- function(y, name) {
  if (length(y) < 2) {
    stop("length of <y> must be at least 2!")
  }
  m <- mean(y)
  std <- stats::sd(y)
  if (std == 0) {
    stop("the varible measurements have zero variance!")
  }
  fun <- function(x) {
    (x - m) / std
  }
  fun_inv <- function(x) {
    x * std + m
  }
  new("lgpscaling", fun = fun, fun_inv = fun_inv, var_name = name)
}

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
        " rows. Found discrepancy between ",
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
