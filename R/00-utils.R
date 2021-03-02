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

#' @rdname likelihood_encoding
is_bin_or_bb <- function(likelihood) {
  likelihood %in% c("binomial", "bb")
}

#' @rdname likelihood_encoding
is_pois_or_nb <- function(likelihood) {
  likelihood %in% c("poisson", "nb")
}

#' Link functions and their inverses
#'
#' @param x the input
#' @param likelihood name of the likelihood model
#' @param a a vector which should be divided
#' elementwise by the vector of numbers of trials
#' @param fit an object of class \linkS4class{lgpfit}
#' @returns transformed input
#' @name link
NULL

#' @rdname link
link <- function(x, likelihood) {
  allowed <- likelihood_list()
  check_allowed(likelihood, allowed)
  if (is_pois_or_nb(likelihood)) {
    x <- log(x)
  } else if (is_bin_or_bb(likelihood)) {
    x <- log(x) - log(1 - x)
  }
  return(x)
}

#' @rdname link
link_inv <- function(x, likelihood) {
  allowed <- likelihood_list()
  check_allowed(likelihood, allowed)
  if (is_pois_or_nb(likelihood)) {
    x <- exp(x)
  } else if (is_bin_or_bb(likelihood)) {
    x <- 1 / (1 + exp(-x))
  }
  return(x)
}

#' @rdname link
divide_by_num_trials <- function(a, fit) {
  check_type(fit, "lgpfit")
  likelihood <- get_obs_model(fit)
  if (!is_bin_or_bb(likelihood)) {
    return(a)
  }
  y_num_trials <- get_num_trials(fit)
  check_lengths(a, y_num_trials)
  a / y_num_trials
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

#' Plot colors to use
#'
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

#' A color palettes
#'
#' @param n an integer from 1 to 5
#' @return an array of \code{n} hex values
#' @family color utilities
#' @name color_palette

#' @rdname color_palette
color_palette <- function(n) {
  c1 <- colorset("brightblue", "mid_highlight")
  c2 <- colorset("red", "mid_highlight")
  c3 <- colorset("orange", "mid_highlight")
  c4 <- colorset("green", "mid_highlight")
  c5 <- colorset("gray", "dark_highlight")
  palette <- c(c1, c2, c3, c4, c5)
  palette[1:n]
}

#' @rdname color_palette
fill_palette <- function(n) {
  c1 <- colorset("brightblue", "mid")
  c2 <- colorset("red", "mid")
  c3 <- colorset("orange", "mid")
  c4 <- colorset("green", "mid")
  c5 <- colorset("gray", "dark")
  palette <- c(c1, c2, c3, c4, c5)
  palette[1:n]
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

#' A color scale.
#'
#' @inheritParams color_palette
#' @return a \code{ggplot} object (if \code{n <= 5}) or NULL (if \code{n > 5})
#' @family color utilities
#' @name scale_color

#' @rdname scale_color
scale_color <- function(n) {
  if (n > 5) {
    return(NULL)
  }
  values <- color_palette(n)
  ggplot2::scale_color_manual(values = values)
}

#' @rdname scale_color
scale_fill <- function(n) {
  if (n > 5) {
    return(NULL)
  }
  values <- fill_palette(n)
  ggplot2::scale_fill_manual(values = values)
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

#' Display a runtime estimation message if starting to analyse a large data set
#'
#' @param num_obs number of observations
#' @param threshold threshold for number of observations
#' @return nothing
large_data_msg <- function(num_obs, threshold) {
  msg <- paste0(
    "WARNING: Number of observations is >= ", threshold,
    ", so sampling can take a long time. See the",
    " 'gradient computation took X seconds' information show by",
    " Stan to estimate total runtime."
  )
  msg <- if (num_obs >= threshold) cat(msg)
}

#' Progress bar for iterative functions
#'
#' \itemize{
#'   \item \code{progbar_header} creates header for a progress bar
#'   \item \code{progbar_print} prints part of the bar depending on iteration
#'   index
#' }
#' @name progbar
NULL

#' @rdname progbar
#' @param L length of bar
progbar_header <- function(L) {
  str <- paste0(seq(10, 100, by = 10))
  a <- formatC(str, width = 3)
  str <- paste0("|  ", a, "%")
  top <- paste(formatC(str, width = 3), collapse = "")
  top <- paste0(top, "|")
  barlen <- nchar(top) - 1
  list(
    header = top,
    idx_print = ceiling(seq(1, L, length.out = barlen))
  )
}

#' @rdname progbar
#' @param idx_print indices of iterations to print a bar block
#' @param idx current iteration index
progbar_print <- function(idx, idx_print) {
  N <- length(which(idx_print == idx))
  str <- paste(rep("=", N), collapse = "")
  cat(str)
}


#' Warning message about using a default prior
#'
#' @param desc parameter description
#' @param name parameter name
#' @param model_desc model description
warn_msg_default_prior <- function(desc, name, model_desc) {
  paste0(
    "Using a default prior for ", desc, " (", name, "), in a model",
    " that ", model_desc, ".",
    " This is not recommended. See the 'Basic usage' tutorial",
    " at https://jtimonen.github.io/lgpr-usage/index.html."
  )
}


#' Get number of nonstationary model components
#'
#' @param stan_input a list containing an element named \code{components}
#' @return an integer
get_num_ns <- function(stan_input) {
  comp <- dollar(stan_input, "components")
  num <- sum(comp[, 5] > 0)
  return(num)
}

#' Repeat a vector as a rows of an array
#'
#' @description Throws an error if \code{v} is \code{NULL}.
#' @param v a vector of length \code{m}
#' @param n number of times to repeat
#' @return returns an array of size \code{n} x \code{m}
#' @family array utilities
repvec <- function(v, n) {
  check_not_null(v)
  m <- length(v)
  A <- matrix(rep(v, n), n, m, byrow = TRUE)
  return(as.array(A))
}

#' Ensure that v is 2-dimensional
#'
#' @param v a vector of length \code{m} or an array of size \code{n} x \code{m}
#' @return returns a 2-dimensional array
#' @family array utilities
ensure_2dim <- function(v) {
  check_not_null(v)
  L <- length(dim(v))
  if (L == 2) {
    out <- v
  } else {
    check_dim(v, 0)
    n <- length(v)
    out <- array(v, dim = c(1, n))
  }
  return(out)
}

#' Compute row variances for a 2-dimensional array
#'
#' @param x an array of size \code{n} x \code{m}
#' @return returns a vector
#' @family array utilities
row_vars <- function(x) {
  check_not_null(x)
  check_dim(x, 2)
  apply(x, 1, stats::var)
}

#' Normalize matrix or data frame rows so that they sum to 1
#'
#' @param x an array of size \code{n} x \code{m}
#' @return an array of size \code{n} x \code{m}
#' @family array utilities
normalize_rows <- function(x) {
  s <- rowSums(x)
  x / s
}

#' Select named row of an array
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param name name of the row
#' @return a vector of length \code{m}
#' @family array utilities
select_row <- function(x, name) {
  df <- data.frame(t(x))
  dollar(df, name)
}

#' Squeeze the second dimension of an array
#'
#' @param x an array of shape \code{n1} x \code{n2} x \code{n3}
#' @return an array of shape \code{n1*n2} x \code{n3}
#' @family array utilities
squeeze_second_dim <- function(x) {
  D_in <- dim(x)
  L <- length(D_in)
  if (L != 3) {
    msg <- paste0("Number of dimensions in <x> must be 3! found = ", L)
    stop(msg)
  }
  if (D_in[2] < 1) {
    msg <- paste0("2nd dimension of <x> must be at least 1! found = ", D_in[2])
    stop(msg)
  }
  D_out <- c(D_in[1] * D_in[2], D_in[3])
  out <- array(0, dim = D_out)
  for (j in seq_len(D_in[2])) {
    i1 <- (j - 1) * D_in[1] + 1
    i2 <- j * D_in[1]
    out[i1:i2, ] <- x[, j, ]
  }
  return(out)
}

#' Array transformation utility
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param L an integer
#' @return a list of length \code{L}, where each element is an array of
#' shape \code{h} x (\code{m}/\code{L}), where \code{h = nrow(x)}
#' @family array utilities
array_to_arraylist <- function(x, L) {
  m <- dim(x)[2]
  if (m %% L) stop("Second dimension of <x> must be a multiple of <L>!")
  out <- list()
  for (k in seq_len(L)) {
    inds <- seq(k, m, by = L)
    out[[k]] <- ensure_2dim(x[, inds])
  }
  return(out)
}

#' Add a vector to each column of an array
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param v a vector of length \code{n}
#' @return an array of shape \code{n} x \code{m}
#' @family array utilities
add_to_columns <- function(x, v) {
  check_dim(x, 2)
  check_length(v, nrow(x))
  D <- ncol(x)
  for (j in seq_len(D)) x[, j] <- x[, j] + v
  return(x)
}

#' Apply possible reduction to rows of an array
#'
#' @param x an array of shape \code{n} x \code{m}
#' @param reduce a function or \code{NULL}
#' @return an array of shape \code{1} x \code{m} if \code{reduce} is not
#' \code{NULL}, otherwise the original array \code{x}
#' @family array utilities
apply_reduce <- function(x, reduce) {
  if (!is.null(reduce)) {
    check_type(reduce, "function")
    x <- apply(x, 2, reduce)
    x <- ensure_2dim(x)
  }
  return(x)
}

#' Return NULL if vector contains only NaN values
#'
#' @param x a vector
#' @return \code{x} unchanged or \code{NULL}
#' @family array utilities
null_if_all_nan <- function(x) {
  L <- length(x)
  S <- sum(is.nan(x))
  out <- if (S < L) x else NULL
  return(out)
}

#' Check that each row of array is identical
#'
#' @param rows rows
#' @return the indentical row
#' @family array utilities
reduce_rows <- function(rows) {
  nam <- rownames(rows)
  R <- dim(rows)[1]
  r1 <- as.numeric(rows[1, ])
  nam1 <- nam[1]
  n <- length(r1)
  if (R == 1) {
    return(r1)
  }
  for (j in 2:R) {
    rj <- as.numeric(rows[j, ])
    namj <- nam[j]
    s <- sum(rj == r1)
    if (s != n) {
      msg <- paste0(
        "For each term with uncrt() or heter() expressions, ",
        "NaNs of the continuous covariate must be on the same",
        " rows. Found discrepancy between ",
        nam1, " and ", namj, "."
      )
      stop(msg)
    }
  }
  return(r1)
}

#' Helper function for creating plot data frames from prediction objects
#'
#' @param x an array with \code{num_draws} rows
#' @return a data frame with \code{num_draws} columns
#' @family array utilities
make_draw_df <- function(x) {
  num_draws <- nrow(x)
  x <- data.frame(t(x))
  colnames(x) <- paste("draw", seq_len(num_draws), sep = "_")
  rownames(x) <- NULL
  return(x)
}

#' Format a 3-dimensional array into a list of 2-dimensional arrays
#'
#' @param x an array
#' @return a list of arrays
#' @family array utilities
arr3_to_list <- function(x) {
  out <- list()
  d <- dim(x)[1]
  for (j in seq_len(d)) {
    out[[j]] <- arr3_select(x, j)
  }
  return(out)
}

#' Select an element along the first dimension of an array, and only drop
#' the first dimension
#'
#' @param x an array of shape\code{d1} x \code{d2} x \code{d3}
#' @param idx the index to select
#' @return an array of shape \code{d2} x \code{d3}
#' @family array utilities
arr3_select <- function(x, idx) {
  DIM <- dim(x)
  array(x[idx, , ], dim = DIM[2:3])
}

#' A safer alternative for the dollar operator
#'
#' @description Requires exact match and throws error if \code{var_name} is
#' not in \code{names(object)}.
#' @param object a list or data frame
#' @param var_name name of the variable to access
#' @return Returns \code{object[[var_name, exact = TRUE]]} if variable
#' exists.
#' @family list utilities
dollar <- function(object, var_name) {
  check_not_null(var_name)
  obj_name <- deparse(substitute(object))
  nams <- names(object)
  if (!(var_name %in% nams)) {
    valid <- paste(nams, collapse = ", ")
    msg <- paste0(
      "Variable with name '", var_name,
      "' not found in <", obj_name, ">!"
    )
    msg <- paste0(msg, " Found: {", valid, "}")
    stop(msg)
  }
  object[[var_name, exact = TRUE]]
}

#' Wrap list into a list of length 1 if the original list is named
#'
#' @param x a list
#' @return a list with no names
#' @family list utilities
list_if_named <- function(x) {
  check_type(x, "list")
  is_named <- !is.null(names(x))
  if (is_named) x <- list(x)
  return(x)
}

#' List elements to matrix rows
#'
#' @param x a list of length \code{m} where each element is a vector of
#' length \code{n}
#' @param n length of each vector
#' @return a matrix with shape \code{m} x \code{n}
#' @family list utilities
list_to_matrix <- function(x, n) {
  m <- length(x)
  A <- array(0, dim = c(m, n))
  for (i in seq_len(m)) {
    A[i, ] <- x[[i]]
  }
  as.matrix(A)
}

#' Matrix rows to a list
#'
#' @param x a matrix or array with \code{m} rows and \code{n} columns
#' @return an unnamed list of length \code{m} where each element is a
#' vector of length \code{n}
#' @family list utilities
matrix_to_list <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- x[i, ]
  }
  return(L)
}

#' Return first list element if the list has length is one
#'
#' @param x a list
#' @return the original list or just its first element if its length is one
#' @family list utilities
simplify_list <- function(x) {
  check_type(x, "list")
  L <- length(x)
  if (L == 1) x <- x[[1]]
  return(x)
}
