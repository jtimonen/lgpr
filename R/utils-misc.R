# Print informational message
log_info <- function(msg, verbose = TRUE) {
  if (verbose) message(msg)
}

# Print progress message
log_progress <- function(msg, verbose = TRUE) {
  if (verbose) cat(msg, "\n")
}

# A safer alternative for the dollar operator
#
# @description Requires exact match and throws an informative error if
# \code{var_name} is not in \code{names(object)}.
# @param object a list or data frame
# @param var_name name of the variable to access
# @return Returns \code{object[[var_name, exact = TRUE]]} if variable
# exists.
dollar <- function(object, var_name) {
  check_not_null(object)
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

# Repeat a vector as a rows of an array
#
# @description Throws an error if \code{v} is \code{NULL}.
# @param v a vector of length \code{m}
# @param n number of times to repeat
# @return returns an array of size \code{n} x \code{m}
repvec <- function(v, n) {
  check_not_null(v)
  m <- length(v)
  A <- matrix(rep(v, n), n, m, byrow = TRUE)
  return(as.array(A))
}

# Ensure that v is 2-dimensional
#
# @param v a vector of length \code{m} or an array of size \code{n} x \code{m}
# @return returns a 2-dimensional array
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


# Squeeze the second dimension of an array
#
# @param x an array of shape \code{n1} x \code{n2} x \code{n3}
# @return an array of shape \code{n1*n2} x \code{n3}
# @family array utilities
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

# Array transformation utility
#
# @param x an array of shape \code{n} x \code{m}
# @param L an integer
# @return a list of length \code{L}, where each element is an array of
# shape \code{h} x (\code{m}/\code{L}), where \code{h = nrow(x)}
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

# Ensure that input is model or fit and return a model
object_to_model <- function(object) {
  if (!is(object, "lgpmodel") && !is(object, "lgpfit")) {
    stop("object must be an lgpmodel or lgpfit object!")
  }
  if (inherits(object, "lgpfit")) {
    out <- object@model
  } else {
    out <- object
  }
  return(out)
}



# List elements to matrix rows
#
# @param x a list of length \code{m} where each element is a vector of
# length \code{n}
# @param n length of each vector
# @return a matrix with shape \code{m} x \code{n}
list_to_matrix <- function(x, n) {
  m <- length(x)
  A <- array(0, dim = c(m, n))
  for (i in seq_len(m)) {
    A[i, ] <- x[[i]]
  }
  as.matrix(A)
}

# Matrix rows to a list
#
# @param x a matrix or array with \code{m} rows and \code{n} columns
# @return an unnamed list of length \code{m} where each element is a
# vector of length \code{n}
matrix_to_list <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- x[i, ]
  }
  return(L)
}

# Repeat a data frame vertically
#
# @param df a data frame
# @param times number of times to repeat
# @return a data frame
rep_df <- function(df, times) {
  row_inds <- rep(seq_len(nrow(df)), times)
  df[row_inds, , drop = FALSE]
}

# Data frame and additional information to long format
#
# @param df a data frame in wide format
# @return a data frame where first column is a factor that determines the
# column variable name in the original data and second column contains
# the actual values
to_long_format <- function(df) {
  nam <- colnames(df)
  x <- c()
  J <- ncol(df)
  N <- nrow(df)
  for (j in seq_len(J)) {
    x <- c(x, as.vector(df[, j]))
  }
  fac <- as.factor(rep(nam, each = N))
  out <- data.frame(fac, x)
  return(out)
}

# Ensure vector has expected length (len) or replicate value to such vector
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

# A linter-friendly way to call a function (fun) with one argument (arg)
call_fun <- function(fun, arg) {
  check_type(fun, "function")
  fun(arg)
}

# Paste function (fun) name and argument (s) name enclosed in parentheses
enclose_fun <- function(s, fun) {
  paste0(fun, "(", s, ")")
}

# Remove quotes and whitespace from a string
simplify_str <- function(s) {
  x <- gsub("[[:space:]]", "", s) # remove whitespace
  x <- gsub("[\",\']", "", x) # remove quotes
  return(x)
}

# Display a runtime estimation message if starting to analyse a large data set
large_data_msg <- function(num_obs, threshold) {
  msg <- paste0(
    "WARNING: Number of observations is >= ", threshold,
    ", so sampling can take a long time. See the",
    " 'gradient computation took X seconds' information shown by",
    " Stan to estimate total runtime."
  )
  show_msg <- (num_obs >= threshold)
  log_info(msg, show_msg)
}

# Create a progress bar header
progbar_setup <- function(L) {
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

# Prints par of progress bar depending on iteration (idx)
progbar_print <- function(idx, idx_print) {
  N <- length(which(idx_print == idx))
  str <- paste(rep("=", N), collapse = "")
  cat(str)
}


# Warning message about using a default prior for (name) when not recommended
warn_msg_default_prior <- function(desc, name, model_desc) {
  paste0(
    "Using a default prior for ", desc, " (", name, "), in a model",
    " that ", model_desc, ".",
    " This is not recommended. See the 'Basic usage' tutorial",
    " at https://jtimonen.github.io/lgpr-usage/index.html."
  )
}

# Compute row variances for a 2-dimensional array
row_vars <- function(x) {
  check_not_null(x)
  check_dim(x, 2)
  apply(x, 1, stats::var)
}

# Normalize matrix or data frame rows so that they sum to 1
normalize_rows <- function(x) {
  s <- rowSums(x)
  x / s
}

# Add vector (v) to each column of 2d array (x)
add_to_columns <- function(x, v) {
  check_dim(x, 2)
  check_length(v, nrow(x))
  D <- ncol(x)
  for (j in seq_len(D)) x[, j] <- x[, j] + v
  return(x)
}

# Apply possible reduction to rows of an array
#
# @param x an array of shape \code{n} x \code{m}
# @param reduce a function or \code{NULL}
# @return an array of shape \code{1} x \code{m} if \code{reduce} is not
# \code{NULL}, otherwise the original array \code{x}
apply_reduce <- function(x, reduce) {
  if (!is.null(reduce)) {
    check_type(reduce, "function")
    x <- apply(x, 2, reduce)
    x <- ensure_2dim(x)
  }
  return(x)
}

# Return NULL if vector contains only NaN values, else original vector (x)
null_if_all_nan <- function(x) {
  L <- length(x)
  S <- sum(is.nan(x))
  out <- if (S < L) x else NULL
  return(out)
}

# Check that each row of array (rows) is identical, an return the row
# Throw error otherwise
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

# Helper function for creating plot data frames from prediction objects
make_draw_df <- function(x) {
  num_draws <- nrow(x)
  x <- data.frame(t(x))
  colnames(x) <- paste("draw", seq_len(num_draws), sep = "_")
  rownames(x) <- NULL
  return(x)
}

# Format a 3-dimensional array into a list of 2-dimensional arrays
# List length is specified by first dimension of x
arr3_to_list <- function(x, list_names) {
  out <- list()
  d <- dim(x)[1]
  for (j in seq_len(d)) {
    out[[j]] <- arr3_select(x, j)
  }
  names(out) <- list_names
  return(out)
}

# Select an element along the first dimension of an array, and only drop
# the first dimension (always returns 2d array)
arr3_select <- function(x, idx) {
  DIM <- dim(x)
  array(x[idx, , ], dim = DIM[2:3])
}


# Wrap list into a list of length 1 if the original list is named
list_if_named <- function(x) {
  check_type(x, "list")
  is_named <- !is.null(names(x))
  if (is_named) x <- list(x)
  return(x)
}

# Return first list element if the list has length is one. Else return the
# whole list unchanged.
simplify_list <- function(x) {
  check_type(x, "list")
  L <- length(x)
  if (L == 1) x <- x[[1]]
  return(x)
}

# Create grouping factor that can be added to a data frame (x)
create_grouping_factor <- function(x, group_by) {
  n <- nrow(x)
  if (is.na(group_by)) {
    x_grp <- as.factor(rep("no grouping", n))
  } else {
    x_grp <- dollar(x, group_by)
    check_type(x_grp, "factor")
  }
  return(x_grp)
}

# Get type of each data frame column except y_name
data_types <- function(data, y_name, verbose) {
  check_type(data, "data.frame")
  data[[y_name]] <- NULL # don't check response variable, it could be integer
  D <- ncol(data)
  nams <- names(data)
  types <- list()
  for (j in seq_len(D)) {
    var <- data[, j]
    if (inherits(var, what = "factor")) {
      types[[j]] <- "factor"
    } else if (inherits(var, what = "numeric")) {
      types[[j]] <- "numeric"
    } else {
      str <- paste(class(var), collapse = ",")
      msg <- paste0(
        "variable ", nams[j], " is neither numeric nor a factor!",
        " (it belongs to classes {", str, "})"
      )
      stop(msg)
    }
    msg <- paste0("Variable type (", nams[j], "): ", types[[j]])
    log_info(msg, verbose)
  }
  names(types) <- nams
  return(types)
}

#' Quick way to create an example lgpfit, useful for debugging
#'
#' @param formula model formula
#' @param likelihood observation model
#' @param chains number of chains to run
#' @param iter number of iterations to run
#' @param num_indiv number of individuals (data simulation)
#' @param num_timepoints number of time points (data simulation)
#' @param ... additional arguments to \code{\link{lgp}}
#' @return An \linkS4class{lgpfit} object created by fitting
#' the example model.
example_fit <- function(formula = y ~ id + age + age | SEX + age | LOC,
                        likelihood = "gaussian",
                        chains = 1,
                        iter = 30,
                        num_indiv = 6,
                        num_timepoints = 5,
                        ...) {
  num_trials <- if (is_bin_or_bb(likelihood)) 20 else NULL
  dat <- simulate_data(
    N = num_indiv, t_data = seq(12, 96, length.out = num_timepoints),
    covariates = c(2, 2),
    relevances = c(0, 1, 0, 1),
    names = c("SEX", "LOC"),
    noise_type = likelihood,
    N_trials = num_trials
  )
  lgp(formula, dat@data,
    likelihood = likelihood,
    chains = chains,
    iter = iter,
    num_trials = num_trials,
    ...
  )
}


# Get total number of post-warmup draws from stanfit
get_num_postwarmup_draws <- function(stan_fit) {
  nrow(as.data.frame(stan_fit))
}
