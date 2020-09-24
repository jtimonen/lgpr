#' Helper function for generic functions
#'
#' @description Helper function for generic functions that work on
#' both of \linkS4class{lgpmodel} and \linkS4class{lgpfit} class objects.
#' @param object an object of class \linkS4class{lgpmodel} or
#' \linkS4class{lgpfit}
#' @return an object of class \linkS4class{lgpmodel}
object_to_model <- function(object) {
  allowed <- c("lgpmodel", "lgpfit")
  check_type(object, allowed)
  if (class(object) == "lgpfit") {
    out <- object@model
  } else {
    out <- object
  }
  return(out)
}

#' Integer encoding of likelihod function names
#'
#' @description
#' \itemize{
#'   \item \code{likelihood_as_int} converts likelihood name to Stan encoding
#'   \item \code{likelihood_as_str} converts the Stan likelihood encoding
#'   to a string
#'   \item \code{likelihood_list} returns the available likelihood names
#' }
#' @param likelihood a string
#' @param index an integer
#' @name likelihood_encoding

#' @rdname likelihood_encoding
likelihood_list <- function() {
  c("gaussian", "poisson", "nb", "binomial", "bb")
}

#' @rdname likelihood_encoding
likelihood_as_str <- function(index) {
  names <- likelihood_list()
  L <- length(names)
  check_interval(index, 1, L)
  name <- names[index]
  return(name)
}

#' @rdname likelihood_encoding
likelihood_as_int <- function(likelihood) {
  likelihood <- tolower(likelihood)
  allowed <- likelihood_list()
  index <- check_allowed(likelihood, allowed)
  return(index)
}

#' Link functions and their inverses
#'
#' @param x the input
#' @param likelihood name of the likelihood model
#' @returns transformed input
#' @name link
NULL

#' @rdname link
link <- function(x, likelihood) {
  allowed <- likelihood_list()
  check_allowed(likelihood, allowed)
  if (likelihood %in% c("poisson", "nb")) {
    x <- log(x)
  } else if (likelihood %in% c("binomial", "bb")) {
    x <- log(x) - log(1 - x)
  }
  return(x)
}

#' @rdname link
link_inv <- function(x, likelihood) {
  allowed <- likelihood_list()
  check_allowed(likelihood, allowed)
  if (likelihood %in% c("poisson", "nb")) {
    x <- exp(x)
  } else if (likelihood %in% c("binomial", "bb")) {
    x <- 1 / (1 + exp(-x))
  }
  return(x)
}

#' Get lgpr version description
#'
#' @export
#' @return package description
get_pkg_description <- function() {
  lgprLib <- dirname(system.file(package = "lgpr"))
  descr <- suppressWarnings(utils::packageDescription("lgpr",
    lib.loc = lgprLib
  ))
  return(descr)
}

#' Get a stan model of the package
#'
#' @export
#' @param name name of the model
#' @return an object of class stanmodel
get_stan_model <- function(name = "lgp") {
  dollar(stanmodels, name)
}

#' Ensure vector has expected length
#'
#' @param len the expected length
#' @param v original vector or just one value that is replicated
#' @return a vector of length \code{len}
ensure_len <- function(v, len) {
  v_name <- deparse(substitute(v))
  L <- length(v)
  if (L == 1) {
    v <- rep(v, len)
  } else if (L != len) {
    msg <- paste0(
      "length of <", v_name, "> was expected to be 1 or ", len,
      ", but found length ", L
    )
    stop(msg)
  }
  return(v)
}

#' Density and quantile functions of the inverse gamma distribution
#'
#' @description Using the same parametrization as Stan. More info
#' \href{https://mc-stan.org/docs/2_24/functions-reference/inverse-gamma-distribution.html}{here}.
#' @param alpha positive real number
#' @param beta positive real number
#' @param x point where to compute the density
#' @param log is log-scale used?
#' @return density/quantile value
#' @name dinvgamma_stanlike
#' @family functions related to the inverse-gamma distribution

#' @rdname dinvgamma_stanlike
dinvgamma_stanlike <- function(x, alpha, beta, log = FALSE) {
  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  if (beta <= 0) {
    stop("beta must be positive")
  }
  t1 <- alpha * log(beta) - lgamma(alpha)
  t2 <- -1 * (alpha + 1) * log(x)
  t3 <- -beta / x
  log_p <- t1 + t2 + t3
  if (log) {
    return(log_p)
  } else {
    return(exp(log_p))
  }
}

#' @param p quantile (must be between 0 and 1)
#' @rdname dinvgamma_stanlike
qinvgamma_stanlike <- function(p, alpha, beta) {
  check_positive(alpha)
  check_positive(beta)
  check_interval(p, 0, 1)
  r <- stats::qgamma(1 - p, shape = alpha, rate = beta)
  return(1 / r)
}

#' Plot colors to use
#'
#' @export
#' @param main Color name. Must be a valid scheme name for
#' \code{\link[bayesplot]{color_scheme_get}}.
#' @param variant Must be one of {"light", "light_highlight", "mid",
#' "mid_highlight", "dark", "dark_highlight"}.
#' @return A hex value of the color.
#' @family color utilities
colorset <- function(main, variant = "mid") {
  scheme <- bayesplot::color_scheme_get(scheme = main)
  col <- scheme[[variant]]
  if (is.null(col)) {
    stop("Invalid color!")
  }
  return(col)
}

#' A color palette function
#'
#' @param n an integer from 1 to 6
#' @return an array of \code{n} hex values
#' @family color utilities
color_palette <- function(n) {
  c1 <- colorset("brightblue", "mid_highlight")
  c2 <- colorset("red", "mid_highlight")
  c3 <- colorset("orange", "mid")
  c4 <- colorset("green", "mid_highlight")
  c5 <- colorset("gray", "dark")
  palette <- c(c1, c2, c3, c4, c5)
  if (n <= 5) {
    out <- palette[1:n]
  } else if (n <= 6) {
    palette <- unlist(bayesplot::color_scheme_get(scheme = "brightblue"))
    out <- as.character(palette)
  } else {
    stop("number of colors can be at most 6")
  }
  return(out)
}

#' Visualize a color palette
#'
#' @inheritParams color_palette
#' @return a \code{ggplot} object
#' @family color utilities
plot_color_palette <- function(n) {
  colors <- color_palette(n)
  x <- rep(c(0, 1), n)
  y <- rep(c(1:n), each = 2)
  col <- as.factor(rep(colors, each = 2))
  df <- data.frame(x, y, col)
  aes <- ggplot2::aes_string(x = x, y = y, color = col, group = col)
  h <- ggplot2::ggplot(df) +
    ggplot2::geom_line(aes, lwd = 1)
  h <- h + ggplot2::scale_color_manual(values = colors)
  blank <- ggplot2::element_blank()
  h <- h + ggplot2::theme(
    axis.text = blank,
    axis.title = blank,
    axis.ticks = blank
  )
  h <- h + ggplot2::theme(legend.position = "none")
  h <- h + ggplot2::ggtitle("Colors")
  return(h)
}

#' A linter-friendly way to call a function
#'
#' @description This exists just to allow calling a function which is a slot
#' of an S4 object, without some linter warnings.
#' @param fun a function which takes one argument
#' @param arg the argument
#' @return The value of \code{fun(arg)}
#' @family function utilities
call_fun <- function(fun, arg) {
  check_type(fun, "function")
  fun(arg)
}

#' Paste function and argument enclosed in parentheses
#'
#' @param s argument name, a string
#' @param fun function name, a string
#' @return a string
#' @family function utilities
enclose_fun <- function(s, fun) {
  paste0(fun, "(", s, ")")
}

#' Remove quotes and whitespace from a string
#'
#' @param s a string
#' @return a string
simplify_str <- function(s) {
  x <- gsub("[[:space:]]", "", s) # remove whitespace
  x <- gsub("[\",\']", "", x) # remove quotes
  return(x)
}
