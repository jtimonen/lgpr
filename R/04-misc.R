#' Character representations of different S4 objects
#'
#' @param x an object of some S4 class
#' @return a character representation of the object
#' @name as_character
NULL

#' @rdname as_character
setMethod("as.character", "lgpexpr", function(x) {
  paste0(x@fun, "(", x@covariate, ")")
})

#' @rdname as_character
setMethod("as.character", "lgpterm", function(x) {
  facs <- x@factors
  desc <- sapply(facs, as.character)
  desc <- paste(desc, collapse = "*")
  return(desc)
})

#' @rdname as_character
setMethod("as.character", "lgpformula", function(x) {
  return(x@call)
})

#' Operations on formula terms and expressions
#'
#' @description
#' \itemize{
#'   \item Sum two \linkS4class{lgprhs}'s to yield an \linkS4class{lgprhs}
#'   \item Sum two \linkS4class{lgpterm}'s to yield an \linkS4class{lgprhs}
#'   \item Sum an \linkS4class{lgprhs} and an \linkS4class{lgpterm}
#'   to yield an \linkS4class{lgprhs}
#'   \item Multiply two \linkS4class{lgpterm}'s to yield
#'   an \linkS4class{lgpterm}
#' }
#' @param e1 The first sum, term or expression
#' @param e2 The second sum, term or expression
#' @name operations
NULL

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgprhs", e2 = "lgprhs"),
  function(e1, e2) {
    new("lgprhs", summands = c(e1@summands, e2@summands))
  }
)

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgprhs", summands = list(e1, e2))
  }
)

#' @rdname operations
setMethod(
  "+", signature(e1 = "lgprhs", e2 = "lgpterm"),
  function(e1, e2) {
    e1 + new("lgprhs", summands = list(e2))
  }
)

#' @rdname operations
setMethod(
  "*", signature(e1 = "lgpterm", e2 = "lgpterm"),
  function(e1, e2) {
    new("lgpterm", factors = c(e1@factors, e2@factors))
  }
)

#' Printing S4 object info using the show generic
#'
#' @name show
#' @param object an object of some S4 class
#' @return the object invisibly
NULL

#' @rdname show
setMethod("show", "lgpformula", function(object) {
  cat(as.character(object))
  invisible(object)
})

#' @rdname show
setMethod("show", "lgprhs", function(object) {
  s <- object@summands
  L <- length(s)
  cat("An object of class lgprhs.\n\n")
  for (j in seq_len(L)) {
    t <- as.character(s[[j]])
    r <- paste0("Term ", j, ": ", t, "\n")
    cat(r)
  }
  invisible(object)
})

#' @rdname show
setMethod("show", "lgpterm", function(object) {
  cat("lgpterm: ")
  cat(as.character(object))
  cat("\n")
  invisible(object)
})

#' @rdname show
setMethod("show", "lgpsim", function(object) {
  msg <- class_info("lgpsim")
  cat(msg)
  invisible(object)
})

#' @rdname show
setMethod("show", "lgpmodel", function(object) {
  msg <- class_info("lgpmodel")
  cat(msg)
  cat("\n")
  model_summary(object)
})

#' @rdname show
setMethod("show", "lgpfit", function(object) {
  msg <- class_info("lgpfit")
  cat(msg)
  cat("\n")
  fit_summary(object)
})

#' @rdname show
setMethod("show", "GaussianPrediction", function(object) {
  fm <- object@f_comp_mean
  num_comps <- length(fm)
  D <- dim(fm[[1]])
  desc <- class_info_pred("GaussianPrediction", num_comps, D)
  cat(desc)
})

#' @rdname show
setMethod("show", "Prediction", function(object) {
  fm <- object@f_comp
  num_comps <- length(fm)
  D <- dim(fm[[1]])
  desc <- class_info_pred("Prediction", num_comps, D)
  cat(desc)
})


#' Class info for show methods
#'
#' @param class_name class name
#' @return a string
class_info <- function(class_name) {
  str <- paste0(
    "An object of class ", class_name, ". See ?",
    class_name, " for more info."
  )
  return(str)
}

#' Class info for show methods of pred objects
#'
#' @param class_name class name
#' @param num_comps number of components
#' @param D dimensions
#' @return a string
class_info_pred <- function(class_name, num_comps, D) {
  desc <- paste0("An object of S4 class '", class_name, "'. ")
  desc <- paste0(desc, "\n - ", num_comps, " components")
  desc <- paste0(desc, "\n - ", D[1], " parameter sets")
  desc <- paste0(desc, "\n - ", D[2], " prediction points")
  return(desc)
}

#' Visualize input warping function with several steepness parameter values
#'
#' @param wrp a vector of values of the warping steepness parameter
#' @param x a vector of input values
#' @param color line color
#' @param alpha line alpha
#' @return a \code{ggplot} object
plot_inputwarp <- function(wrp,
                           x,
                           color = colorset("red", "dark"),
                           alpha = 0.5) {
  x <- sort(x)
  L <- length(x)
  S <- length(wrp)
  W <- matrix(0, S, L)
  for (i in seq_len(S)) {
    w <- cpp_warp_input(x, a = wrp[i])
    W[i, ] <- w
  }
  af <- as.factor(rep(1:S, each = L))
  df <- data.frame(rep(x, S), as.vector(t(W)), af)
  colnames(df) <- c("x", "w", "idx")

  # Create ggplot object
  aes <- ggplot2::aes_string(x = "x", y = "w", group = "idx")
  plt <- ggplot2::ggplot(df, aes)

  # Add titles and labels
  plt <- plt + ggplot2::labs(x = "Input", y = "Warped input") +
    ggplot2::ggtitle("Input-warping function")
  plt <- plt + ggplot2::ylim(-1.0, 1.0)

  # Plot the actual lines of interest
  plt <- plt + ggplot2::geom_line(color = color, alpha = alpha)
  return(plt)
}
