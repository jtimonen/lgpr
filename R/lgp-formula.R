#' Create a model formula
#'
#' @param formula an object of class \code{formula}
#' @return an object of class \code{\link{lgpformula}}
lgp_formula <- function(formula) {
  c_str <- class(formula)
  if (c_str != "formula") {
    stop("<formula> must have class formula! found = ", c_str)
  }
  f_str <- as.character(formula)
  if (length(f_str) != 3) {
    stop("Invalid formula: as.character(formula) should have length 3!")
  }
  text <- f_str[3]
  components <- eval(parse(text = text)) # TODO: allow a + a:b syntax
  components <- ensure_lgpsum(components)

  new("lgpformula",
    response = f_str[2],
    components = components,
    call = f_str
  )
}

#' Promote lgpterm or lgpproduct to lgpsum
#'
#' @param x an object of class \code{\link{lgpterm}}, \code{\link{lgpproduct}}
#' or \code{\link{lgpsum}}
#' @return an object of class \code{\link{lgpterm}}
ensure_lgpsum <- function(x) {
  if (class(x) == "lgpterm") {
    x <- new("lgpproduct", factors = list(x))
  }
  if (class(x) == "lgpproduct") {
    x <- new("lgpsum", summands = list(x))
  }
  return(x)
}

#' Create a stationary GP term
#'
#' @param covariate name of a continuous covariate
#' @return an object of class \code{\link{lgpterm}}
gp <- function(covariate) {
  covariate <- noquote(deparse1(substitute(covariate)))
  new("lgpterm",
    covariate = as.character(covariate),
    type = "continuous",
    fun = "gp"
  )
}

#' Create a nonstationary GP term
#'
#' @param covariate name of a continuous covariate
#' @return an object of class \code{\link{lgpterm}}
gp_ns <- function(covariate) {
  covariate <- noquote(deparse1(substitute(covariate)))
  new("lgpterm",
    covariate = as.character(covariate),
    type = "continuous",
    fun = "gp_ns"
  )
}

#' Create a zerosum term
#'
#' @param covariate name of a categorical covariate
#' @return an object of class \code{\link{lgpterm}}
zerosum <- function(covariate) {
  covariate <- noquote(deparse1(substitute(covariate)))
  new("lgpterm",
    covariate = as.character(covariate),
    type = "discrete",
    fun = "zerosum"
  )
}

#' Create a categorical term
#'
#' @param covariate name of a categorical covariate
#' @return an object of class \code{\link{lgpterm}}
categorical <- function(covariate) {
  covariate <- noquote(deparse1(substitute(covariate)))
  new("lgpterm",
    covariate = as.character(covariate),
    type = "discrete",
    fun = "categorical"
  )
}

#' Create a mask term
#'
#' @param covariate name of a categorical covariate
#' @return an object of class \code{\link{lgpterm}}
mask <- function(covariate) {
  covariate <- noquote(deparse1(substitute(covariate)))
  new("lgpterm",
    covariate = as.character(covariate),
    type = "discrete",
    fun = "mask"
  )
}

# Creating a product term
setMethod(
  "*", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgpproduct", factors = list(e1, e2))
  }
)

# Summing two products
setMethod(
  "+", signature(e1 = "lgpproduct", e2 = "lgpproduct"),
  function(e1, e2) {
    new("lgpsum", summands = list(e1, e2))
  }
)

# Summing two terms
setMethod(
  "+", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    c1 <- new("lgpproduct", factors = list(e1))
    c2 <- new("lgpproduct", factors = list(e2))
    new("lgpsum", summands = list(c1, c2))
  }
)

# Summing product and term
setMethod(
  "+", signature(e1 = "lgpproduct", e2 = "lgpterm"),
  function(e1, e2) {
    c <- new("lgpproduct", factors = list(e2))
    new("lgpsum", summands = list(e1, c))
  }
)

# Summing term and product
setMethod(
  "+", signature(e1 = "lgpterm", e2 = "lgpproduct"),
  function(e1, e2) {
    c <- new("lgpproduct", factors = list(e1))
    new("lgpsum", summands = list(c, e2))
  }
)

# Summing sum and product
setMethod(
  "+", signature(e1 = "lgpsum", e2 = "lgpproduct"),
  function(e1, e2) {
    a <- e1@summands
    a[[length(a) + 1]] <- e2
    new("lgpsum", summands = a)
  }
)

# Summing sum and term
setMethod(
  "+", signature(e1 = "lgpsum", e2 = "lgpterm"),
  function(e1, e2) {
    c <- new("lgpproduct", factors = list(e2))
    e1 + c
  }
)
