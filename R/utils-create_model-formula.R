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
#' @inheritParams create_model.likelihood
#' @return an object of class \linkS4class{lgpformula}
#' @family internal model creation functions
create_model.formula <- function(formula, data, verbose = FALSE) {
  log_progress("Parsing formula...", verbose)
  advanced <- is_advanced_formula(formula)
  if (!advanced) formula <- formula_to_advanced(formula, data, verbose)
  fp <- as.character(formula)
  formula_str <- paste(fp[2], fp[1], fp[3])
  log_info(paste0("Formula interpreted as: ", formula_str), verbose)
  parse_formula_advanced(formula)
}

# Check if formula is in advanced format
is_advanced_formula <- function(formula) {
  f <- split_formula(formula)
  rhs <- dollar(f, "right")
  out <- grepl("(", rhs, fixed = TRUE)
  return(out)
}

# Split formula into left and right-hand side
split_formula <- function(formula) {
  check_type(formula, "formula")
  f_str <- as.character(formula)
  check_length(f_str, 3)
  out <- list(left = f_str[2], right = f_str[3])
  return(out)
}

# Translate formula to advanced syntax format
formula_to_advanced <- function(formula, data, verbose) {
  f <- split_formula(formula)
  rhs <- dollar(f, "right")
  y_name <- dollar(f, "left")
  rhs <- rhs_to_advanced(rhs, data, verbose, y_name)
  f_str <- paste(y_name, "~", rhs)
  out <- stats::formula(f_str)
  return(out)
}

# Translate formula right-hand side
rhs_to_advanced <- function(rhs, data, verbose, y_name) {
  types <- data_types(data, y_name, verbose)
  terms <- rhs_to_terms(rhs)
  rhs_adv <- ""
  idx <- 0
  for (t in terms) {
    idx <- idx + 1
    covs <- strsplit(t, split = "|", fixed = TRUE)[[1]]
    if (length(covs) == 1) {
      check_in_data(covs[1], data, "data")
      t1 <- dollar(types, covs[1])
      if (t1 == "numeric") {
        f1 <- "gp"
      } else if (t1 == "factor") {
        f1 <- "zs"
      } else {
        stop("Invalid data type at this point. Please report a bug.")
      }
      term <- enclose_fun(covs[1], f1)
    } else {
      check_length(covs, 2)
      check_in_data(covs[1], data, "data")
      check_in_data(covs[2], data, "data")
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

# Split a formula right-hand side to terms
rhs_to_terms <- function(rhs) {
  x <- simplify_str(rhs)
  out <- strsplit(x, split = "+", fixed = TRUE)[[1]]
  return(out)
}

# Parse a model formula that uses the advanced formula syntax and
# return an object of class \linkS4class{lgpformula}
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

# Parse a string representation of the right-hand side of an advanced formula
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

# Parse a string representation of one advanced formula term
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

# Parse a string representation of one advanced formula expression
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

# Get names of all variables appearing in a term of advanced formula
term_variables <- function(term) {
  a <- character()
  for (f in term@factors) {
    a <- c(a, f@covariate)
  }
  return(a)
}

# Get names of all variables appearing on advanced formula right-hand side
rhs_variables <- function(rhs) {
  a <- character()
  for (s in rhs@summands) {
    a <- c(a, term_variables(s))
  }
  return(a)
}
