#' Create a model formula
#'
#' @param formula The model formula, where
#' \itemize{
#'   \item it must contain exatly one tilde (\code{~}), with response
#'   variable on the left-hand side and model terms on the right-hand side
#'   \item terms are be separated by a plus (\code{+}) sign
#'   \item all variables appearing in \code{formula} must be
#'   found in \code{data}
#' }
#' See the ´Syntax for formula terms´ section in the documentation of
#' \code{\link{lgp}} for instructions on how to specify the model terms.
#' @return an object of class \linkS4class{lgpformula}
#' @family model formula functions
parse_formula <- function(formula) {

  # Check input
  check_type(formula, "formula")
  f_str <- as.character(formula)
  check_length(f_str, 3)
  rhs_str <- f_str[3]

  # Check if formula is in advanced format and translate if not
  advanced <- is_advanced_formula(rhs_str)
  if (!advanced) rhs_str <- formula_to_advanced(rhs_str)

  # Parse the right-hand side and create lgpformula object
  terms <- parse_rhs(rhs_str)
  new("lgpformula",
    y_name = f_str[2],
    terms = terms,
    call = paste(f_str[2], f_str[1], f_str[3])
  )
}

#' Parse a model formula that uses the common formula syntax
#'
#' @param description Translates formula to advanced syntax format.
#' @inheritParams is_advanced_formula
#' @return an object of class \linkS4class{lgpformula}
#' @family model formula functions
formula_to_advanced <- function(rhs) {
  terms <- rhs_to_terms(rhs)
  for (t in terms) {
    print(t)
  }
}

#' Check if formula is in advanced format.
#'
#' @param rhs The formula right-hand side in text format.
#' @return TRUE or FALSE
#' @family model formula functions
is_advanced_formula <- function(rhs) {
  grepl("(", rhs, fixed = TRUE)
}

#' Split a formula to terms
#'
#' @inheritParams is_advanced_formula
#' @return a character vector
#' @family model formula functions
rhs_to_terms <- function(rhs) {
  x <- parse_whitespace(rhs)
  strsplit(x, split = "+", fixed = TRUE)[[1]]
}

#' Parse a string representation of the right-hand side of a formula
#'
#' @inheritParams is_advanced_formula
#' @return an object of class \linkS4class{lgprhs}
#' @family model formula functions
parse_rhs <- function(rhs) {
  terms <- rhs_to_terms(rhs)
  D <- length(terms)
  f <- parse_term(terms[1])
  if (D == 1) {
    out <- new("lgprhs", summands = list(f))
    return(out)
  } else {
    for (j in 2:D) {
      f <- f + parse_term(terms[j])
    }
    return(f)
  }
}

#' Parse a string representation of one formula term
#'
#' @param term a term without any plus signs, e.g. \code{b*c}
#' @return an object of class \linkS4class{lgpterm}
#' @family model formula functions
parse_term <- function(term) {
  factors <- strsplit(term, split = "*", fixed = TRUE)[[1]] # split to factors
  D <- length(factors)
  if (D > 4) {
    stop("One formula term can have at most four expressions! found = ", D)
  }
  e <- parse_expr(factors[[1]])
  f <- new("lgpterm", factors = list(e))
  if (D == 1) {
    return(f)
  } else {
    for (j in 2:D) {
      e <- parse_expr(factors[[j]])
      f <- f * new("lgpterm", factors = list(e))
    }
    return(f)
  }
}

#' Parse a string representation of one formula expression
#'
#' @param expr the expression as a string, e.g. \code{"gp(x)"}
#' @return an object of class \linkS4class{lgpexpr}
#' @family model formula functions
parse_expr <- function(expr) {
  num_open <- lengths(regmatches(expr, gregexpr("[(]", expr)))
  num_close <- lengths(regmatches(expr, gregexpr("[)]", expr)))
  if (num_open == 0) {
    msg <- paste0(
      "Zero opening brackets found in the expression <",
      expr, ">. Be sure not to mix the advanced and simple",
      " formula syntaxes. See ?lgp."
    )
    stop(msg)
  }
  ok_paren <- (num_open == 1) && (num_close == 1)
  if (!ok_paren) {
    msg <- paste0(
      "Each expression must contain exactly one opening and ",
      "closing parenthesis! Found ", num_open, " opening and ",
      num_close, " closing brackets in the expression <", expr, ">."
    )
    stop(msg)
  }
  parsed <- strsplit(expr, "[()]")[[1]]
  check_length(parsed, 2)
  out <- new("lgpexpr", covariate = parsed[2], fun = parsed[1])
  return(out)
}

#' Remove quotes and whitespace from a formula
#'
#' @inheritParams is_advanced_formula
#' @return a character string
#' @family model formula functions
parse_whitespace <- function(rhs) {
  x <- gsub("[[:space:]]", "", rhs) # remove whitespace
  x <- gsub("[\",\']", "", x) # remove quotes
  return(x)
}

#' Get names of all variables appearing in a term
#'
#' @param term an object of class \linkS4class{lgpterm}
#' @return a list of variable names
#' @family model formula functions
term_variables <- function(term) {
  a <- character()
  for (f in term@factors) {
    a <- c(a, f@covariate)
  }
  return(a)
}

#' Get names of all variables appearing on formula right-hand side
#'
#' @param rhs an object of class \linkS4class{lgprhs}
#' @return a list of variable names
#' @family model formula functions
rhs_variables <- function(rhs) {
  a <- character()
  for (s in rhs@summands) {
    a <- c(a, term_variables(s))
  }
  return(a)
}
