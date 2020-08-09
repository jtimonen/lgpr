# Summing two rhs's
setMethod(
  "+", signature(e1 = "lgprhs", e2 = "lgprhs"),
  function(e1, e2) {
    new("lgprhs", summands = c(e1@summands, e2@summands))
  }
)

# Summing two terms
setMethod(
  "+", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgprhs", summands = list(e1, e2))
  }
)

# Summing rhs and term
setMethod(
  "+", signature(e1 = "lgprhs", e2 = "lgpterm"),
  function(e1, e2) {
    e1 + new("lgprhs", summands = list(e2))
  }
)

# Multiplying of two exprs
setMethod(
  "*", signature(e1 = "lgpexpr", e2 = "lgpexpr"),
  function(e1, e2) {
    new("lgpterm", factors = list(e1, e2))
  }
)

#' Parse a string representation of one formula expression
#'
#' @param expr the expression as a string, e.g. \code{"gp(x)"}
#' @return an object of class \code{\link{lgpexpr}}
parse_expr <- function(expr) {
  num_open <- lengths(regmatches(expr, gregexpr("[(]", expr)))
  num_close <- lengths(regmatches(expr, gregexpr("[)]", expr)))
  if (num_open != 1) {
    stop(
      "Each expression must contain exactly one opening parenthesis! ",
      "Found ", num_open, " in <", expr, ">."
    )
  }
  if (num_close != 1) {
    stop(
      "Each expression must contain exactly one closing parenthesis! ",
      "Found ", num_close, " in <", expr, ">."
    )
  }
  parsed <- strsplit(expr, "[()]")[[1]]
  L <- length(parsed)
  if (L > 2) {
    stop("invalid expression ", expr)
  }
  out <- new("lgpexpr", covariate = parsed[2], fun = parsed[1])
  return(out)
}

#' Parse a string representation of one formula term
#'
#' @param term a term without any plus signs, e.g. \code{b*c}
#' @return an object of class \code{\link{lgpterm}}
parse_term <- function(term) {
  factors <- strsplit(term, split = "*", fixed = TRUE)[[1]] # split to factors
  D <- length(factors)
  if (D > 2) {
    stop("One formula term can have at most two factors! found = ", D)
  }
  f1 <- parse_expr(factors[1])
  if (D == 2) {
    f2 <- parse_expr(factors[2])
    return(f1 * f2)
  } else {
    out <- new("lgpterm", factors = list(f1))
    return(out)
  }
}

#' Parse a string representation of the right-hand side of a formula
#'
#' @param rhs the formula right-hand side, e.g. \code{a + b*c}
#' @return an object of class \code{\link{lgprhs}}
parse_rhs <- function(rhs) {
  x <- gsub("[[:space:]]", "", rhs) # remove whitespace
  x <- gsub("[\",\']", "", x) # remove quotes
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

#' Create a model formula
#'
#' @param formula A formula specified using the common \code{\link{formula}}
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
#' @return an object of class \code{\link{lgpformula}}
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
  new("lgpformula",
    response = f_str[2],
    terms = terms,
    call = paste(f_str[2], f_str[1], f_str[3])
  )
}

#' Get names of all variables appearing in a term
#'
#' @param term an object of class \code{lgpterm}
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
#' @param rhs an object of class \code{lgprhs}
#' @return a list of variable names
rhs_variables <- function(rhs) {
  a <- character()
  for (s in rhs@summands) {
    a <- c(a, term_variables(s))
  }
  return(a)
}
