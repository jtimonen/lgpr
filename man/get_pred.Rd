% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main-get_pred.R
\name{get_pred}
\alias{get_pred}
\title{Extract model predictions and function posteriors}
\usage{
get_pred(fit, draws = NULL, reduce = NULL, verbose = TRUE)
}
\arguments{
\item{fit}{An object of class \linkS4class{lgpfit}.}

\item{draws}{Indices of parameter draws to use, or \code{NULL} to use all
draws.}

\item{reduce}{Reduction for parameters draws. Can be a function that
is applied to reduce all parameter draws into one parameter set, or
\code{NULL} (no reduction). Has no effect if \code{draws} is specified.}

\item{verbose}{Should more information and a possible progress bar be
printed?}
}
\value{
an object of class \linkS4class{GaussianPrediction} or
\linkS4class{Prediction}
}
\description{
\emph{NOTE:} It is not recommended for users to call this. Use
\code{\link{pred}} instead.
}
