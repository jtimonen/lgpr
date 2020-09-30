#' Create a model formula
#'
#' @description Checks if formula is in advanced format and translates if not.
#' @param formula The model formula, where
#' \itemize{
#'   \item it must contain exatly one tilde (\code{~}), with response
#'   variable on the left-hand side and model terms on the right-hand side
#'   \item terms are be separated by a plus (\code{+}) sign
#'   \item all variables appearing in \code{formula} must be
#'   found in \code{data}
#' }
#' See the "Model formula syntax" section below (\code{\link{lgp}}) for
#' instructions on how to specify the model terms.
#' @inheritParams parse_response
#' @param verbose Should more verbose output be printed?
#' @return an object of class \linkS4class{lgpformula}
#' @family model formula functions
parse_formula <- function(formula, data, verbose = FALSE) {
  advanced <- is_advanced_formula(formula)
  if (!advanced) formula <- formula_to_advanced(formula, data)
  if (verbose) cat("Formula converted to:\n  ")
  if (verbose) print(formula)
  parse_formula_advanced(formula)
}

#' Check if formula is in advanced format.
#'
#' @description Checks if right-hand side contains opening parentheses.
#' @inheritParams parse_formula
#' @return TRUE or FALSE
#' @family model formula functions
is_advanced_formula <- function(formula) {
  f <- split_formula(formula)
  rhs <- dollar(f, "right")
  out <- grepl("(", rhs, fixed = TRUE)
  return(out)
}

#' Split formula into left and right-hand side
#'
#' @inheritParams parse_formula
split_formula <- function(formula) {
  check_type(formula, "formula")
  f_str <- as.character(formula)
  check_length(f_str, 3)
  out <- list(left = f_str[2], right = f_str[3])
  return(out)
}

#' Parse a model formula that uses the common formula syntax
#'
#' @description Translates formula to advanced syntax format.
#' @inheritParams is_advanced_formula
#' @inheritParams parse_formula
#' @return an object of class \linkS4class{lgpformula}
#' @name formula_to_advanced
NULL

#' @rdname formula_to_advanced
formula_to_advanced <- function(formula, data) {
  f <- split_formula(formula)
  rhs <- dollar(f, "right")
  y_name <- dollar(f, "left")
  rhs <- rhs_to_advanced(rhs, data)
  f_str <- paste(y_name, "~", rhs)
  out <- stats::formula(f_str)
  return(out)
}

#' @rdname formula_to_advanced
#' @inheritParams rhs_to_terms
rhs_to_advanced <- function(rhs, data) {
  types <- data_types(data)
  terms <- rhs_to_terms(rhs)
  rhs_adv <- ""
  idx <- 0
  for (t in terms) {
    idx <- idx + 1
    covs <- strsplit(t, split = "|", fixed = TRUE)[[1]]
    if (length(covs) == 1) {
      check_in_data(covs[1], data)
      t1 <- dollar(types, covs[1])
      f1 <- if ("factor" %in% t1) "zs" else "gp"
      term <- enclose_fun(covs[1], f1)
    } else {
      check_length(covs, 2)
      check_in_data(covs[1], data)
      check_in_data(covs[2], data)
      e1 <- enclose_fun(covs[1], "gp")
      t2 <- dollar(types, covs[2])
      f2 <- if ("factor" %in% t2) "zs" else "gp"
      e2 <- enclose_fun(covs[2], f2)
      term <- paste0(e1, "*", e2)
    }
    rhs_adv <- if (idx == 1) term else paste(rhs_adv, "+", term)
  }
  return(rhs_adv)
}


#' Split a formula to terms
#'
#' @param rhs the formula right-hand side in text format
#' @return a character vector
#' @family model formula functions
rhs_to_terms <- function(rhs) {
  x <- simplify_str(rhs)
  out <- strsplit(x, split = "+", fixed = TRUE)[[1]]
  return(out)
}
