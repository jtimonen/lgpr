% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lgp.R
\name{create_model}
\alias{create_model}
\title{Create a model}
\usage{
create_model(
  formula,
  data,
  likelihood = "gaussian",
  prior = NULL,
  c_hat = NULL,
  num_trials = NULL,
  options = NULL,
  verbose = FALSE
)
}
\arguments{
\item{formula}{A formula specified using the common \code{formula}
syntax, such as \code{y ~ x1 + x2:z1 + } \code{x2:z2 + z2}.
\itemize{
  \item The formula must contain exatly one tilde (\code{~}), with response
  variable on the left-hand side and model terms on the right-hand side.
  \item Terms are be separated by a plus (\code{+}) sign.
  \item Terms can consits of a single variable name, such as \code{x}, or
  an interaction of two variables, such as \code{x:z}.
  \item In single-variable terms, the variable can be either continuous or
  categorical, whereas in interaction terms the variable
  on the left-hand side of the colon (\code{:}) has to be continuous and the
  one on the right-hand side has to be categorical (a factor).
  \item All variables appearing in \code{formula} must be
  found in \code{data}.
}}

\item{data}{A \code{data.frame} where each column corresponds to one
variable, and each row is one observation. Continuous covariates and the
response variable must have type \code{"numeric"} and categorical covariates
must have type \code{"factor"}. Missing values should be indicated with
\code{NaN} or \code{NA}. The response variable cannot contain missing
values.}

\item{likelihood}{Determines the observation model. Must be either
\code{"gaussian"} (default), \code{"poisson"}, \code{"nb"} (negative
binomial) or \code{"binomial"}. To use Bernoulli likelihood, use
\code{likelihood="binomial"} and set \code{num_trials} as a vector of ones.}

\item{prior}{A named list, defining the prior distribution of model
(hyper)parameters.}

\item{c_hat}{The GP mean. Must be a vector of length \code{dim(data)[1]}, or
a real number defining a constant GP mean. If \code{NULL}, this is set to
 \itemize{
   \item \code{c_hat = 0}, if \code{likelihood} is \code{"gaussian"}, because
   with Gaussian likelihood the response variable is by default centered to
   have zero mean.
   \item \code{c_hat = } \code{log(mean(y))} if \code{likelihood} is
   \code{"poisson"} or \code{"nb"},
   \item \code{c_hat = } \code{log(p/(1-p))}, where
   \code{p = mean(y/num_trials)} if \code{likelihood} is \code{"binomial"},
 }
where \code{y} denotes the response variable. You can modify this vector to
account for normalization between data points. With Gaussian likelihood
though, do not modify this argument, normalize the data beforehand instead.}

\item{num_trials}{This argument (number of trials) is only needed when
likelihood is binomial. Must have length one or equal to number of data
points. Setting \code{num_trials=1} corresponds to Bernoulli observation
model.}

\item{options}{A named list with the following possible fields:
\itemize{
  \item \code{sample_f} Determines if the function values are be sampled
  (must be \code{TRUE} if likelihood is not \code{"gaussian"}).
  \item \code{skip_generated} If this is true, the generated quantities
  block of Stan is skipped.
  \item \code{delta} Amount of added jitter to ensure positive definite
  covariance matrices.
}}

\item{verbose}{Should more verbose output be printed?}
}
\description{
Create a model
}