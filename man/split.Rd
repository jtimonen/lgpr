% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main-datautil.R
\name{split}
\alias{split}
\alias{split_by_factor}
\alias{split_within_factor}
\alias{split_within_factor_random}
\alias{split_random}
\alias{split_data}
\title{Split data into training and test sets}
\usage{
split_by_factor(data, test, var_name = "id")

split_within_factor(data, idx_test, var_name = "id")

split_within_factor_random(data, k_test = 1, var_name = "id")

split_random(data, p_test = 0.2, n_test = NULL)

split_data(data, i_test, sort_ids = TRUE)
}
\arguments{
\item{data}{a data frame}

\item{test}{the levels of the factor that will be used as test data}

\item{var_name}{name of a factor in the data}

\item{idx_test}{indices point indices with the factor}

\item{k_test}{desired number of test data points per each level of the
factor}

\item{p_test}{desired proportion of test data}

\item{n_test}{desired number of test data points (if NULL, \code{p_test}
is used to compute this)}

\item{i_test}{test data row indices}

\item{sort_ids}{should the test indices be sorted into increasing order}
}
\value{
a named list with names \code{train}, \code{test}, \code{i_train}
and \code{i_test}
}
\description{
\itemize{
  \item \code{split_by_factor} splits according to given factor
  \item \code{split_within_factor} splits according to given
  data point indices within the same level of a factor
  \item \code{split_within_factor_random} selects k points
  from each level of a factor uniformly at random as test data
  \item \code{split_random} splits uniformly at random
  \item \code{split_data} splits according to given data rows
}
}
\seealso{
Other data frame handling functions: 
\code{\link{add_dis_age}()},
\code{\link{add_factor}()},
\code{\link{add_factor_crossing}()},
\code{\link{adjusted_c_hat}()},
\code{\link{new_x}()}
}
\concept{data frame handling functions}
