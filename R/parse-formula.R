#' Create a model formula
#'
#' @param formula A formula specified using the common \code{formula}
#' syntax, such as \code{y ~ x1 + x2:z1 + } \code{x2:z2 + z2}.
#' \itemize{
#'   \item The formula must contain exatly one tilde (\code{~}), with response
#'   variable on the left-hand side and model terms on the right-hand side.
#'   \item Terms are be separated by a plus (\code{+}) sign.
#'   \item Terms can consits of a single variable name, such as \code{x}, or
#'   an interaction of two variables, such as \code{x:z}.
#'   \item In single-variable terms, the variable can be either continuous or
#'   categorical, whereas in interaction terms the variable
#'   on the left-hand side of the colon (\code{:}) has to be continuous and the
#'   one on the right-hand side has to be categorical (a factor).
#'   \item All variables appearing in \code{formula} must be
#'   found in \code{data}.
#' }
#' @return an object of class \linkS4class{lgpformula}
parse_formula <- function(formula) {
  c_str <- class(formula)
  if (c_str != "formula") {
    stop("<formula> must have class formula! found = ", c_str)
  }
  f_str <- as.character(formula)
  if (length(f_str) != 3) {
    stop("Invalid formula: as.character(formula) should have length 3!")
  }
  text <- f_str[3]
  terms <- parse_rhs(text)
  call_str <- paste(f_str[2], f_str[1], f_str[3])
  new("lgpformula",
    y_name = f_str[2],
    terms = terms,
    call = parse_whitespace(call_str)
  )
}

#' Parse a string representation of the right-hand side of a formula
#'
#' @inheritParams parse_whitespace
#' @return an object of class \linkS4class{lgprhs}
parse_rhs <- function(rhs) {
  x <- parse_whitespace(rhs)
  terms <- strsplit(x, split = "+", fixed = TRUE)[[1]] # split to terms
  D <- length(terms)
  if (D < 1) {
    stop("<rhs> must have at least one term!")
  }
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
parse_expr <- function(expr) {
  num_open <- lengths(regmatches(expr, gregexpr("[(]", expr)))
  num_close <- lengths(regmatches(expr, gregexpr("[)]", expr)))
  ok_paren <- (num_open == 1) && (num_close == 1)
  if (!ok_paren) {
    msg <- paste0(
      "Each expression must contain exactly one opening and ",
      "closing parenthesis! Found {", num_open, ", ",
      num_close, "} in <", expr, ">."
    )
    stop(msg)
  }
  parsed <- strsplit(expr, "[()]")[[1]]
  L <- length(parsed)
  if (L > 2) {
    stop("invalid expression ", expr)
  }
  out <- new("lgpexpr", covariate = parsed[2], fun = parsed[1])
  return(out)
}


#' Remove quotes and whitespace from a formula
#'
#' @param rhs the formula right-hand side, e.g. \code{a + b*c}
#' @return a character string
parse_whitespace <- function(rhs) {
  x <- gsub("[[:space:]]", "", rhs) # remove whitespace
  x <- gsub("[\",\']", "", x) # remove quotes
  return(x)
}





#' Get names of all variables appearing in a term
#'
#' @param term an object of class \linkS4class{lgpterm}
#' @return a list of variable names
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
rhs_variables <- function(rhs) {
  a <- character()
  for (s in rhs@summands) {
    a <- c(a, term_variables(s))
  }
  return(a)
}
