
#' Validate S4 class objects
#'
#' @param object an object to validate
#' @return \code{TRUE} if valid, otherwise reasons for invalidity
#' @name validate
NULL

#' @rdname validate
check_lgpexpr <- function(object) {
  errors <- character()
  v0 <- nchar(object@covariate) > 0
  valid_funs <- c(
    "gp", "gp_warp", "gp_warp_vm", "categ", "zerosum", "heter", "uncrt"
  )
  v1 <- object@fun %in% valid_funs
  if (!v0) {
    errors <- c(errors, "covariate name cannot be empty")
  }
  if (!v1) {
    str <- paste0(valid_funs, collapse = ", ")
    msg <- paste0(
      "<fun> must be one of {", str,
      "}, found = '", object@fun, "'"
    )
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
}

#' @rdname validate
check_lgpformula <- function(object) {
  covs <- rhs_variables(object@terms)
  r <- object@y_name
  errors <- character()
  if (r %in% covs) {
    msg <- "the response variable cannot be also a covariate"
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
}

#' @rdname validate
check_lgpscaling <- function(object) {
  a <- 1.2321
  a_mapped <- object@fun(a)
  b <- object@fun_inv(a_mapped)
  diff <- abs(a - b)
  errors <- character()
  if (diff > 1e-6) {
    err <- paste("<f_inv> is not an inverse function of <f>, diff = ", diff)
    errors <- c(errors, err)
  }
  if (nchar(object@var_name) < 1) {
    err <- "variable name length must be at least 1 character"
    errors <- c(errors, err)
  }
  out <- if (length(errors) > 0) errors else TRUE
  return(out)
}
