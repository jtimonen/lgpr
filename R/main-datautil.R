
#' Set the GP mean vector, taking TMM or other normalization
#' into account
#'
#' @export
#' @description Creates the \code{c_hat} input for \code{lgp},
#' so that it accounts for normalization between data points in the
#' \code{"poisson"} or \code{"nb"} observation model
#' @param y response variable, vector of length \code{n}
#' @param norm_factors normalization factors, vector of length \code{n}
#' @return a vector of length \code{n}, which can be used as
#' the \code{c_hat} input to the \code{lgp} function
#' @family data frame handling functions
adjusted_c_hat <- function(y, norm_factors) {
  check_lengths(y, norm_factors)
  check_non_negative_all(y)
  check_integer_all(y)
  check_positive_all(norm_factors)
  c_hat <- log(mean(y))
  c_hat <- rep(c_hat, length(y))
  c_hat <- c_hat + log(norm_factors)
  return(c_hat)
}

#' Easily add a categorical covariate to a data frame
#'
#' @export
#' @param data the original data frame
#' @param x A named vector containing the category for each individual.
#' The names should specify the individual id.
#' @param id_var name of the id variable in \code{data}
#' @return A data frame with one column added. The new column will
#' have same name as the variable passed as input \code{x}.
#' @family data frame handling functions
add_factor <- function(data, x, id_var = "id") {
  data <- convert_to_data_frame(data)
  name <- deparse(substitute(x))
  bad <- name %in% colnames(data)
  if (bad) stop("<data> already contains a variable called '", name, "'!")
  check_named(x)
  x_id <- as.numeric(names(x))
  data_id <- dollar(data, id_var)
  uid <- unique(x_id)
  new_factor <- rep(0, length(data_id))
  for (id in uid) {
    i_data <- which(data_id == id)
    i_new <- which(x_id == id)
    new_factor[i_data] <- x[i_new]
  }
  data[[name]] <- as.factor(new_factor)
  return(data)
}

#' Add a crossing of two factors to a data frame
#'
#' @param data a data frame
#' @param fac1 name of first factor, must be found in \code{df}
#' @param fac2 name of second factor, must be found in \code{df}
#' @param new_name name of the new factor
#' @return a data frame
#' @family data frame handling functions
add_factor_crossing <- function(data, fac1, fac2, new_name) {
  df <- convert_to_data_frame(data)
  a <- dollar(df, fac1)
  b <- dollar(df, fac2)
  check_not_null(new_name)
  check_type(a, "factor")
  check_type(b, "factor")
  df[[new_name]] <- interaction(a, b, sep = "*")
  return(df)
}

#' Easily add the disease-related age variable to a data frame
#'
#' @export
#' @description Creates the disease-related age covariate vector based on the
#' disease initiation times and adds it to the data frame
#' @param data the original data frame
#' @param t_init A named vector containing the observed initiation or onset
#' time for each individual. The names, i.e. \code{names(t_init)}, should
#' specify the individual id.
#' @param id_var name of the id variable in \code{data}
#' @param time_var name of the time variable in \code{data}
#' @return A data frame with one column added. The new column will
#' be called \code{dis_age}. For controls, its value will be \code{NaN}.
#' @family data frame handling functions
add_dis_age <- function(data, t_init, id_var = "id", time_var = "age") {
  data <- convert_to_data_frame(data)
  bad <- "dis_age" %in% colnames(data)
  if (bad) stop("<data> already contains a variable called 'dis_age'!")
  check_named(t_init)
  x_id <- as.numeric(names(t_init))
  data_id <- dollar(data, id_var)
  data_age <- dollar(data, time_var)
  uid <- unique(x_id)
  dage <- rep(NaN, length(data_id))
  for (id in uid) {
    i_data <- which(data_id == id)
    i_new <- which(x_id == id)
    dage[i_data] <- data_age[i_data] - t_init[i_new]
  }
  data$dis_age <- dage
  return(data)
}

#' Split data into training and test sets
#'
#' @name split
#' @param data a data frame
#' @param var_name name of a factor in the data
#' @description
#' \itemize{
#'   \item \code{split_by_factor} splits according to given factor
#'   \item \code{split_within_factor} splits according to given
#'   data point indices within the same level of a factor
#'   \item \code{split_within_factor_random} selects k points
#'   from each level of a factor uniformly at random as test data
#'   \item \code{split_random} splits uniformly at random
#'   \item \code{split_data} splits according to given data rows
#' }
#' @return a named list with names \code{train}, \code{test}, \code{i_train}
#' and \code{i_test}
#' @family data frame handling functions
NULL

#' @rdname split
#' @param test the levels of the factor that will be used as test data
split_by_factor <- function(data, test, var_name = "id") {
  check_type(data, "data.frame")
  fac <- dollar(data, var_name)
  check_type(fac, "factor")
  i_test <- which(fac %in% test)
  split_data(data, i_test)
}

#' @rdname split
#' @param idx_test indices point indices with the factor
split_within_factor <- function(data, idx_test, var_name = "id") {
  data <- convert_to_data_frame(data)
  id <- dollar(data, var_name)
  check_type(id, "factor")
  uid <- unique(id)
  i_test <- c()
  for (i in uid) {
    inds <- which(id == i)
    i_test <- c(i_test, inds[idx_test])
  }
  split_data(data, i_test)
}

#' @rdname split
#' @param k_test desired number of test data points per each level of the
#' factor
split_within_factor_random <- function(data, k_test = 1, var_name = "id") {
  data <- convert_to_data_frame(data)
  id <- dollar(data, var_name)
  check_type(id, "factor")
  uid <- unique(id)
  i_test <- c()
  for (i in uid) {
    inds <- which(id == i)
    L <- length(inds)
    i_sel <- sample.int(L, k_test)
    i_test <- c(i_test, inds[i_sel])
  }
  split_data(data, i_test)
}

#' @rdname split
#' @param p_test desired proportion of test data
#' @param n_test desired number of test data points (if NULL, \code{p_test}
#' is used to compute this)
split_random <- function(data, p_test = 0.2, n_test = NULL) {
  data <- convert_to_data_frame(data)
  n_total <- dim(data)[1]
  if (is.null(n_test)) {
    n_test <- round(p_test * n_total)
  }
  i_test <- sample.int(n_total, size = n_test)
  split_data(data, i_test)
}

#' @rdname split
#' @param i_test test data row indices
#' @param sort_ids should the test indices be sorted into increasing order
split_data <- function(data, i_test, sort_ids = TRUE) {
  data <- convert_to_data_frame(data)
  i_test <- if (sort_ids) sort(i_test, decreasing = FALSE) else sort_ids
  n_total <- dim(data)[1]
  i_train <- setdiff(seq_len(n_total), i_test)

  # Return
  list(
    train = data[i_train, ],
    test = data[i_test, ],
    i_train = i_train,
    i_test = i_test
  )
}

#' Create prediction points
#'
#' @description Replaces a continuous variable \code{x} in the data frame, and
#' possibly another continuous variable \code{x_ns} derived from it, with new
#' values, for each level of a grouping factor (usually id)
#' @export
#' @param data A data frame. Can also be an \linkS4class{lgpfit} or
#' \linkS4class{lgpmodel} object, in which case data is extracted from it.
#' @param group_by name of the grouping variable, must be a factor
#' in \code{data} (or use \code{group_by=NA} to create a dummy grouping
#' factor which has only one value)
#' @param x of the variable along which to extend,
#' must be a numeric in \code{data}
#' @param x_ns of a nonstationary variable derived from \code{x},
#' must be a numeric in \code{data}
#' @param x_values the values of \code{x} to set for each individual
#' @return a data frame containing the following columns
#' \itemize{
#'  \item all factors in the original \code{data}
#'  \item \code{x}
#'  \item \code{x_ns} (unless it is NULL)
#' }
#'
#' @family data frame handling functions
new_x <- function(data, x_values, group_by = "id", x = "age", x_ns = NULL) {
  data <- allow_data_model_or_fit(data)
  data <- convert_to_data_frame(data)
  check_not_null(x)
  check_not_null(x_values)
  check_in_data(x, data, "data")
  if (is.na(group_by)) {
    x_grp <- create_grouping_factor(data, group_by) # util
    data["group__"] <- x_grp
    group_by <- "group__"
  } else {
    check_in_data(group_by, data, "data")
  }
  df <- pick_one_row_each(data, group_by)
  k <- length(x_values)
  col_names <- if (is.null(x_ns)) x else c(x, x_ns)
  df <- select_factors_and(df, col_names)
  N <- nrow(df)
  inds <- rep(seq_len(N), each = k)
  df <- df[inds, ]
  x_val_rep <- rep(x_values, times = N)
  df[[x]] <- x_val_rep
  if (!is.null(x_ns)) {
    t0 <- get_teff_obs(data, group_by, x, x_ns)
    t0 <- as.numeric(t0)
    df[[x_ns]] <- x_val_rep - rep(t0, each = k)
  }
  rownames(df) <- NULL
  return(df)
}

# Get observed effect times from a data frame
get_teff_obs <- function(data, group_by = "id", x = "age",
                         x_ns = "diseaseAge") {
  check_type(data, "data.frame")
  df <- pick_one_row_each(data, group_by)
  times <- dollar(df, x) - dollar(df, x_ns)
  names(times) <- dollar(df, group_by)
  return(times)
}

# For each unique value of a factor fac, pick one row from data
pick_one_row_each <- function(data, fac) {
  check_type(data, "data.frame")
  z <- dollar(data, fac)
  check_type(z, "factor")
  rows <- c()
  for (lev in levels(z)) {
    inds <- which(z == lev)
    rows <- c(rows, inds[1])
  }
  data[rows, ]
}

# Select data columns which are factors or have name specified by valid
select_factors_and <- function(data, valid) {
  check_type(data, "data.frame")
  col_inds <- c()
  D <- ncol(data)
  nams <- names(data)
  for (j in seq_len(D)) {
    a <- is.factor(data[, j])
    b <- nams[j] %in% valid
    if (a || b) {
      col_inds <- c(col_inds, j)
    }
  }
  data[, col_inds]
}

#' Vizualizing longitudinal data
#'
#' @export
#' @param data A data frame.
#' @param x_name Name of x-axis variable.
#' @param y_name Name of the y-axis variable.
#' @param group_by Name of grouping variable (must be a factor).
#' @param color_by Name of coloring variable (must be a factor).
#' @param facet_by Name of the faceting variable (must be a factor).
#' @param highlight Value of category of the \code{group_by} variable
#' that is highlighted. Can only be used if \code{color_by} is \code{NULL}.
#' @param main main plot title
#' @param sub plot subtitle
#' @return a \code{\link[ggplot2]{ggplot}} object
#' @family data frame handling
plot_data <- function(data,
                      x_name = "age",
                      y_name = "y",
                      group_by = "id",
                      facet_by = NULL,
                      color_by = NULL,
                      highlight = NULL,
                      main = NULL,
                      sub = NULL) {
  data <- allow_data_model_or_fit(data)
  data <- convert_to_data_frame(data)

  # Check that needed variables exist
  check_in_data(x_name, data, "data")
  check_in_data(y_name, data, "data")
  check_in_data(group_by, data, "data")

  # Create initial plot and add data
  df <- data[c(x_name, y_name, group_by, color_by, facet_by)]
  df <- plot_data_add_highlight_factor(df, group_by, highlight)
  skip_hl <- is.null(highlight)
  color_by <- if (skip_hl) color_by else paste0(group_by, "_")
  aes <- plot_data_aes(x_name, y_name, group_by, color_by)
  h <- ggplot2::ggplot(df, aes)

  # Add data
  h <- h + ggplot2::geom_line() + ggplot2::geom_point()

  # Add titles, faceting and coloring
  titles <- plot_data_titles(main, sub, data, group_by)
  label <- dollar(titles, "main")
  subtitle <- dollar(titles, "sub")
  h <- h + ggplot2::ggtitle(label = label, subtitle = subtitle)
  if (!is.null(facet_by)) {
    f <- stats::as.formula(paste("~", facet_by))
    h <- h + ggplot2::facet_wrap(f, labeller = ggplot2::label_both)
  }

  return(h)
}

# Create aes for plot_data
plot_data_aes <- function(x_name, y_name, group_by, color_by) {
  if (is.null(color_by)) {
    aes <- ggplot2::aes_string(x = x_name, y = y_name, group = group_by)
  } else {
    aes <- ggplot2::aes_string(
      x = x_name, y = y_name, group = group_by,
      color = color_by
    )
  }
  return(aes)
}

# Create titles for plot_data
plot_data_titles <- function(main, sub, data, group_by) {
  g <- data[[group_by]]
  N <- length(unique(g))
  n <- dim(data)[1]
  if (is.null(main)) {
    main <- paste0(n, " data points")
  }
  if (is.null(sub)) {
    sub <- paste0(
      "Points that share the same value for '", group_by,
      "' are connected by a line (", N, " levels)."
    )
  }
  list(main = main, sub = sub, N = N, n = n)
}

# Add factor to data frame for highlighting in plot
plot_data_add_highlight_factor <- function(df, group_by, highlight) {
  if (!is.null(highlight)) {
    check_length(highlight, 1)
    g <- df[[group_by]]
    s <- sum(g == highlight)
    if (s == 0) {
      str <- paste(levels(g), collapse = ", ")
      msg <- paste0(
        "Invalid <highlight> argument ", highlight, "! The ",
        "possible values of ", group_by, " are: {", str, "}."
      )
      stop(msg)
    }
    hl <- 1 + as.numeric(g == highlight)
    levels <- c("other", highlight)
    name <- paste0(group_by, "_")
    df[[name]] <- as.factor(levels[hl])
  }
  return(df)
}

# This allows data utilities to work with lgpmodel or lgpfit, too
allow_data_model_or_fit <- function(data) {
  if (is(data, "lgpfit")) data <- get_data(get_model(data))
  if (is(data, "lgpmodel")) data <- get_data(data)
  check_type(data, "data.frame")
  return(data)
}

# Convert data to data.frame
convert_to_data_frame <- function(data) {
  check_type(data, "data.frame") # check if inherits from data.frame
  cname <- class(data)
  if (length(cname) > 1) {
    warning("data is not a plain data.frame, converting using as.data.frame()")
    data <- as.data.frame(data)
  }
  data
}



#' Function for reading the built-in proteomics data
#'
#' @export
#' @param parentDir Path to local parent directory for the data.
#' If this is \code{NULL}, data is downloaded from
#' \url{https://github.com/jtimonen/lgpr-usage/tree/master/data/proteomics}.
#' @param protein Index or name of protein.
#' @param verbose Can this print some output?
#' @return a \code{data.frame}
#' @family built-in data
read_proteomics_data <- function(parentDir = NULL, protein = NULL,
                                 verbose = TRUE) {
  if (is.null(protein)) stop("specify name or index of protein!")
  REPO_PATH <- "jtimonen/lgpr-usage/master/data/proteomics/"
  RAW_PATH <- paste0("https://raw.githubusercontent.com/", REPO_PATH)
  fn_X <- "liu_preproc_X.csv"
  fn_Y <- "liu_preproc_Y.csv"
  if (is.null(parentDir)) {
    msg <- "Given parentDir was NULL, downloading the data from internet..."
    log_info(msg, verbose)
    fn_X <- paste0(RAW_PATH, fn_X)
    fn_Y <- paste0(RAW_PATH, fn_Y)
    X_data <- read_csv_url(fn_X, verbose, header = TRUE, sep = ",")
    Y_data <- read_csv_url(fn_Y, verbose, header = TRUE, sep = ",")
  } else {
    fn_X <- paste0(parentDir, "/", fn_X)
    fn_Y <- paste0(parentDir, "/", fn_Y)
    X_data <- utils::read.csv(fn_X, header = TRUE, sep = ",")
    Y_data <- utils::read.csv(fn_Y, header = TRUE, sep = ",")
  }

  # Get protein
  names <- colnames(Y_data)
  if (!is.character(protein)) {
    pname <- names[protein]
  } else {
    pname <- protein
  }
  msg <- paste0("Read data for protein '", pname, "'.")
  log_info(msg, verbose)

  # Remove rows containing NaNs for the protein
  y <- Y_data[[pname]]
  notnan <- which(!is.nan(y))
  num_nan <- length(which(is.nan(y)))
  data <- data.frame(cbind(X_data, y))
  data <- data[notnan, ]
  msg <- paste0(
    "Removed ", num_nan, " rows with NaN value as the ",
    "protein measurement."
  )
  if (num_nan > 0) log_info(msg, verbose)

  # Categorical variables to factors
  data$id <- as.factor(data$id)
  data$group <- as.factor(data$group)
  data$sex <- as.factor(data$sex)
  levels(data$group)[levels(data$group) == 0] <- "Control"
  levels(data$group)[levels(data$group) == 1] <- "Case"
  levels(data$sex)[levels(data$sex) == 0] <- "Female"
  levels(data$sex)[levels(data$sex) == 1] <- "Male"
  return(data)
}

# Read a csv from URL
read_csv_url <- function(fn, verbose, ...) {
  msg <- paste0("Reading ", fn, " with ?accessType=DOWNLOAD")
  log_info(msg, verbose)
  url <- paste0(fn, "?accessType=DOWNLOAD")
  text <- RCurl::getURL(url)
  utils::read.csv(text = text, ...)
}
