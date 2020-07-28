
#' Get values of function components at data points
#' @param df A \code{stanfit} object as data frame,
#' obtained as \code{as.data.frame(stanfit)}
#' @param model An object of class \code{lgpmodel}
#' @return An array of size
#' \code{n_samples} x \code{n_data} x \code{n_components+2}
get_function_components_from_df_all <- function(df, model) {
  n_samples <- dim(df)[1]
  n <- model@stan_dat$n
  d <- sum(model@stan_dat$D)
  res <- array(0, c(n_samples, n, d + 2))
  for (i in 1:n_samples) {
    res[i, , ] <- get_function_components_from_df(df[i, ], model)
  }
  return(res)
}


#' Get values of function components at data points, for one MCMC sample
#' @param pars A data frame representing one parameter sample,
#' i.e one row of \code{as.data.frame(stanfit)}, where stanfit
#' is an object of class \code{stanfit}
#' @param model An object of class \code{lgpmodel}
#' @return A matrix of size \code{n_data} x \code{n_components+2}
get_function_components_from_df <- function(pars, model) {
  F_sampled <- as.numeric(model@stan_dat$F_IS_SAMPLED)
  n <- model@stan_dat$n
  d <- sum(model@stan_dat$D)
  pn <- names(pars)
  if (!F_sampled) {
    i1 <- which(grepl("F_mean_cmp", pn))
    i2 <- which(grepl("F_mean_tot", pn))
    F_cmp <- matrix(as.numeric(pars[i1]), n, d, byrow = TRUE)
    F_tot <- matrix(as.numeric(pars[i2]), n, 1, byrow = TRUE)
  } else {
    i1 <- which(grepl("F", pn))
    F_cmp <- matrix(as.numeric(pars[i1]), n, d, byrow = TRUE)
    F_tot <- matrix(rowSums(F_cmp), n, 1, byrow = TRUE)
  }
  G_tot <- get_g_from_f(F_tot, model)
  res <- cbind(F_cmp, F_tot, G_tot)
  colnames(res) <- c(model@info$component_names, "f", "g")
  return(res)
}

#' Get a posterior estimate of model (hyper)parameters
#'
#' @param object An (incomplete) object of class \code{lgpfit}.
#' @param type Must be "mean", "median", or "map".
#' @return a data frame
hyperparam_estimate <- function(object, type = "mean") {
  DF <- as.data.frame(object@stan_fit)
  nam <- names(DF)
  i_eta <- which(grepl("ETA", nam))
  i_f <- which(grepl("F", nam))
  i_lp <- which(grepl("lp__", nam))
  if (type == "mean") {
    DF <- colMeans(DF)
  } else if (type == "median") {
    DF <- apply(DF, 2, stats::median)
  } else if (type == "map") {
    i_max <- which(DF$lp__ == max(DF$lp__))
    DF <- colMeans(DF[i_max, ]) # colmeans just to get correct data type
  } else {
    stop("Invalid type!")
  }
  DF <- DF[-c(i_f, i_eta, i_lp)]
  return(DF)
}


#' Get a set of model (hyper)parameter samples
#'
#' @param object An (incomplete) object of class \code{lgpfit}.
#' @param samples Sample indices. If NULL, all samples are taken.
#' @return a data frame
hyperparam_samples <- function(object, samples = NULL) {
  nam <- names(object@stan_fit)
  i_eta <- which(grepl("ETA", nam))
  i_f <- which(grepl("F", nam))
  i_lp <- which(grepl("lp__", nam))
  nam <- nam[-c(i_f, i_eta, i_lp)]
  ext <- rstan::extract(object@stan_fit, pars = nam)
  n <- length(ext[[1]])
  cn <- names(ext)
  d <- length(cn)
  OUT <- matrix(0, n, d)
  for (j in 1:d) {
    OUT[, j] <- ext[[j]]
  }
  OUT <- data.frame(OUT)
  colnames(OUT) <- cn
  if (is.null(samples)) {
    samples <- 1:n
  }
  OUT <- OUT[samples, ]
  return(OUT)
}


#' Get average runtime of a chain
#' @param object An object of class \code{lgpfit}.
#' @return Average runtimes for warmup and sampling
get_runtime <- function(object) {
  TIM <- rstan::get_elapsed_time(object@stan_fit)
  n_chains <- dim(TIM)[1]
  t1 <- round(mean(TIM[, 1]), 2)
  t2 <- round(mean(TIM[, 2]), 2)
  return(list(warmup = t1, sampling = t2))
}

#' Get signal on data scale from process f
#' @param f A vector
#' @param model an object of class \code{lgpmodel}
#' @return A vector g
get_g_from_f <- function(f, model) {
  sdat <- model@stan_dat
  LH <- sdat$LH
  if (LH == 1 || LH == 0) {
    g <- f
  } else if (LH == 2 || LH == 3) {
    g <- exp(f + sdat$C_hat)
  } else if (LH == 4) {
    logistic <- function(x) {
      1 / (1 + exp(-x))
    }
    g <- sdat$N_trials * logistic(sdat$C_hat + f)
  } else {
    stop("Unknown likelihood!")
  }
  return(g)
}
