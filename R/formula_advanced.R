#' Parse a model formula that uses the advanced formula syntax
#'
#' @inheritParams parse_formula
#' @return an object of class \linkS4class{lgpformula}
#' @family advanced formula parsers
parse_formula_advanced <- function(formula) {
  f <- split_formula(formula)
  rhs <- dollar(f, "right")
  terms <- parse_rhs(rhs)
  y_name <- dollar(f, "left")
  new("lgpformula",
    y_name = y_name,
    terms = terms,
    call = paste(y_name, "~", rhs)
  )
}

#' Parse a string representation of the right-hand side of a formula
#'
#' @inheritParams rhs_to_terms
#' @return an object of class \linkS4class{lgprhs}
#' @family advanced formula parsers
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
#' @family advanced formula parsers
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
#' @family advanced formula parsers
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

#' Get names of all variables appearing in a term
#'
#' @param term an object of class \linkS4class{lgpterm}
#' @return a list of variable names
#' @family advanced formula parsers
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
#' @family advanced formula parsers
rhs_variables <- function(rhs) {
  a <- character()
  for (s in rhs@summands) {
    a <- c(a, term_variables(s))
  }
  return(a)
}
