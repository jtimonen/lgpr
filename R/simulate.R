#' Generate an artificial longitudinal data set
#'
#' @export
#' @description Generate an artificial longitudinal data set.
#' @inheritParams sim_create_x
#' @inheritParams sim_create_f
#' @inheritParams sim_create_y
#' @param N_affected Number of diseased individuals that are affected by the
#' disease. This defaults to the number of diseased individuals. This argument
#' can only be given if \code{covariates} contains a zero.
#' @param t_observed Determines how the disease effect time is observed. This
#' can be any function that takes the real disease effect time as an argument
#' and returns the (possibly randomly generated) observed onset/initiation time.
#' Alternatively, this can be a string of the form \code{"after_n"} or
#' \code{"random_p"} or \code{"exact"}.
#' @param f_var variance of f
#' @param c_hat A constant added to f
#' @param verbose Verbosity mode.
#' @return An object of class \linkS4class{lgpsim}.
#' @examples
#' # Generate Gaussian data
#' dat <- simulate_data(N = 4, t_data = c(6, 12, 24, 36, 48), snr = 3)
#'
#' # Generate negative binomially distributed count data
#' dat <- simulate_data(
#'   N = 6, t_data = seq(2, 10, by = 2), noise_type = "nb",
#'   phi = 2
#' )
simulate_data <- function(N,
                          t_data,
                          covariates = c(),
                          names = NULL,
                          relevances = c(1, 1, rep(1, length(covariates))),
                          n_categs = rep(2, sum(covariates %in% c(2, 3))),
                          t_jitter = 0,
                          lengthscales =
                            rep(12, 2 + sum(covariates %in% c(0, 1, 2))),
                          f_var = 1,
                          noise_type = "gaussian",
                          snr = 3,
                          phi = 1,
                          N_affected = round(N / 2),
                          t_effect_range = "auto",
                          t_observed = "after_0",
                          c_hat = 0,
                          dis_fun = "gp_vm",
                          bin_kernel = FALSE,
                          steepness = 0.5,
                          vm_params = c(0.025, 1),
                          continuous_info = list(
                            mu = c(pi / 8, pi, -0.5),
                            lambda = c(pi / 8, pi, 1)
                          ),
                          N_trials = 1,
                          verbose = FALSE,
                          force_zeromean = TRUE) {

  # Input checks
  noise_type <- tolower(noise_type)
  if (N < 2) stop("There must be at least 2 individuals!")
  if (length(t_data) < 3) {
    stop("There must be at least 3 time points per individual!")
  }
  names <- sim_check_covariates(covariates, relevances, names, n_categs)
  if (N_affected > round(N / 2)) {
    stop("N_affected cannot be greater than round(N/2)!")
  }
  if (is.unsorted(covariates)) {
    stop("The covariates vector must be increasing!")
  }

  # Generate the covariates
  IN <- sim_create_x(
    N = N,
    covariates = covariates,
    names = names,
    n_categs = n_categs,
    t_data = t_data,
    t_jitter = t_jitter,
    t_effect_range = t_effect_range,
    continuous_info = continuous_info
  )
  if (verbose) {
    cat(IN$info)
  }

  # Compute X_affected
  k <- length(t_data)
  X_affected <- c(rep(1, N_affected * k), rep(0, (N - N_affected) * k))

  # Simulate the components F
  COMP <- sim_create_f(
    IN$X,
    covariates,
    relevances,
    lengthscales,
    X_affected,
    dis_fun,
    bin_kernel,
    steepness,
    vm_params,
    force_zeromean
  )
  FFF <- COMP$FFF

  # Total signal f is a (scaled) sum of the components plus a constant
  f <- rowSums(FFF)
  SD <- stats::sd(f)
  if (SD > 0) {
    f <- sqrt(f_var) / SD * f
  } else {
    # if the sum of components is zero, define total signal as just noise
    f <- stats::rnorm(n = length(f))
  }
  f <- f + c_hat

  # Generate noise
  NOISY <- sim_create_y(noise_type, f,
    snr = snr, phi = phi,
    N_trials = N_trials
  )
  g <- NOISY$g
  y <- NOISY$y
  noise <- y - g

  # Create the output objects
  dat <- cbind(IN$X, y)
  rownames(dat) <- 1:(N * k)
  comp <- cbind(FFF, f, g, noise, y)

  # Real diseaseAges to observed ones
  OBSERVED <- sim_data_to_observed(dat, t_observed)

  # Signal and noise sums of squares
  SSR <- sum((g - mean(g))^2)
  SSE <- sum(noise^2)

  # Create rest of the fields
  teff_true <- IN$onsets
  teff_obs <- OBSERVED$onsets_observed
  teff <- list(true = teff_true, observed = teff_obs)

  info <- list(
    par_ell = lengthscales,
    par_cont = IN$par_cont,
    p_signal = SSR / (SSR + SSE)
  )

  # Return S4 class object
  new("lgpsim",
    data = OBSERVED$dat,
    response = "y",
    components = comp,
    kernel_matrices = COMP$KKK,
    effect_times = teff,
    info = info
  )
}


#' Generate noisy observations

#' @param noise_type Either "gaussian", "poisson", "nb" (negative binomial) or
#' "binomial".
#' @param snr The desired signal-to-noise ratio. This argument is valid
#' only with \cr
#' \code{noise_type = "gaussian"}.
#' @param phi The dispersion parameter for negative binomial data. The variance
#' is g + g^2/phi.
#' @param N_trials The number of trials parameter for binomial data.
#' @param f The underlying signal.
#' @return A list \code{out}, where
#' \itemize{
#'   \item \code{out$g} is \code{f} mapped through an inverse link function and
#'   \item \code{out$y} is the noisy response variable.
#' }
sim_create_y <- function(noise_type, f, snr, phi, N_trials) {
  L <- length(f)
  y <- rep(0, L)
  if (noise_type == "gaussian") {
    sf <- stats::var(f)
    sigma_n <- 1 / sqrt(snr) * sqrt(sf)
    y_n <- stats::rnorm(n = L, mean = 0, sd = 1)
    y_n <- sigma_n * (y_n - mean(y_n)) / stats::sd(y_n)
    g <- f
    y <- g + y_n
  } else if (noise_type == "poisson") {
    g <- exp(f)
    for (i in 1:L) {
      y[i] <- stats::rpois(n = 1, lambda = g[i])
    }
  } else if (noise_type == "nb") {
    g <- exp(f)
    for (i in 1:L) {
      y[i] <- stats::rnbinom(n = 1, mu = g[i], size = phi)
    }
  } else if (noise_type == "binomial") {
    g <- 1 / (1 + exp(-f))
    for (i in 1:L) {
      y[i] <- stats::rbinom(n = 1, size = N_trials, prob = g[i])
    }
  } else {
    stop(paste("Invalid input noise_type = '", noise_type, "'", sep = ""))
  }
  NOISY <- list(g = g, y = y)
  return(NOISY)
}


#' Input check for the covariates-related arguments of \code{simulate_data}
#' @param covariates argument to \code{simulate_data}
#' @param relevances argument to \code{simulate_data}
#' @param names argument to \code{simulate_data}
#' @param n_cat the \code{n_categs} argument to \code{simulate_data}
#' @return the covariate names
sim_check_covariates <- function(covariates, relevances, names, n_cat) {

  # Validity check
  d0 <- sum(covariates == 0)
  if (d0 > 1) stop("There can be only one diseaseAge component!")
  L1 <- length(covariates)
  L2 <- length(relevances)
  L3 <- sum(covariates %in% c(2, 3))
  if ((L1 + 2) != L2) stop("length(relevances) must be length(covariates) + 2")
  if (L3 != length(n_cat)) stop("The argument n_cat has invalid length!")

  if (!is.null(names)) {
    # Check input names
    L4 <- length(names)
    if (L4 != L1) stop("length(names) must be the same as length(covariates)")
    names <- c("id", "age", names)
  } else {
    # Generate default names
    names <- sim_generate_names(covariates)
  }

  if (sum(relevances < 0) > 0) {
    stop("The relevances must be non-negative!")
  }

  return(names)
}


#' Generate names for covariates
#' @param covariates vector of covariate types
#' @return covariate names
sim_generate_names <- function(covariates) {


  # Get default names
  names <- c("id", "age")
  def <- c("x", "z", "offset", "group")
  D <- sim_create_x_D(covariates)
  if (D[1] == 1) {
    names <- c(names, "diseaseAge")
  }
  if (D[2] > 0) {
    if (D[2] > 1) {
      names <- c(names, paste(def[1], 1:D[2], sep = ""))
    } else {
      names <- c(names, def[1])
    }
  }
  if (D[3] > 0) {
    if (D[3] > 1) {
      names <- c(names, paste(def[2], 1:D[3], sep = ""))
    } else {
      names <- c(names, def[2])
    }
  }
  if (D[4] > 0) {
    if (D[4] > 1) {
      names <- c(names, paste(def[3], 1:D[4], sep = ""))
    } else {
      names <- c(names, def[3])
    }
  }
  if (D[5] > 0) {
    names <- c(names, "group")
  }

  return(names)
}


#' Real generated disease ages to observed ones
#' @param dat data frame
#' @param t_observed Determines how the disease onset is observed. See
#' documentation of \code{\link{simulate_data}}.
#' @return a new data frame and observed onsets
sim_data_to_observed <- function(dat, t_observed) {
  flag <- !("diseaseAge" %in% colnames(dat))
  id <- dat$id
  uid <- unique(id)
  N <- length(uid)
  onsets_observed <- rep(NaN, N)
  names(onsets_observed) <- c(1:N)
  if (flag) {
    # No modifications to data frame needed
    ret <- list(dat = dat, onsets_observed = onsets_observed)
    return(ret)
  } else {
    age <- dat$age
    dis_age <- dat$diseaseAge
    j <- 0
    for (ID in uid) {
      j <- j + 1
      inds <- which(id == ID)
      age_i <- age[inds]
      dag_i <- dis_age[inds]
      if (is.nan(dag_i[1])) {
        # not a diseased individual
      } else {

        # how many points are there after the real onset?
        irem <- which(dag_i > 0)
        rem <- age_i[irem]
        t0_real <- -dag_i[1] + age_i[1]

        # sample the observed onset t0
        if (is.function(t_observed)) {
          t_possible <- t_observed(t0_real)
          inds0 <- which(rem > t_possible)
          if (length(inds0) < 1) {
            stop("There are no data points after t = ", t_possible, "!")
          } else {
            idx0 <- inds0[1] # next meas. point after detection is possible
            t0 <- rem[idx0]
          }
        } else {
          parsed <- sim_parse_t_obs(t_observed)
          if (parsed$type == "random") {
            idx0 <- sim_rtgeom(length(irem), parsed$value)
            if (idx0 > length(rem)) {
              stop("Not enough data points to go that far!")
            }
            t0 <- rem[idx0]
          } else if (parsed$type == "after") {
            idx0 <- parsed$value + 1
            if (idx0 > length(rem)) {
              stop("Not enough data points to go that far!")
            }
            t0 <- rem[idx0]
          } else {
            t0 <- t0_real
          }
        }
        onsets_observed[j] <- t0

        # update the diseaseAge covariate
        dat$diseaseAge[inds] <- age_i - t0
      }
    }
    ret <- list(dat = dat, onsets_observed = onsets_observed)
    return(ret)
  }
}

#' Sample from the 'truncated geometric' distribution
#' @param p a number between 0 and 1
#' @param s an integer
#' @param n number of samples
#' @return an integer from the interval 1...n
sim_rtgeom <- function(s, p, n = 1) {
  s_seq <- 0:(s - 1)
  prob <- p^s_seq
  r <- sample.int(n = s, size = n, prob = prob, replace = TRUE)
  return(r)
}


#' Parse the t_observed argument of \code{simulate_data}
#' @param t_observed a string
#' @return a list with a name and number
sim_parse_t_obs <- function(t_observed) {
  parts <- strsplit(t_observed, "_")[[1]]
  if (parts[1] == "after") {
    type <- "after"
  } else if (parts[1] == "random") {
    type <- "random"
  } else if (parts[1] == "exact") {
    type <- "exact"
  } else {
    stop("unknown keyword ", parts[1])
  }
  value <- as.numeric(parts[2])
  return(list(type = type, value = value))
}

#' Simulate latent function components for longitudinal data analysis
#'
#' @param X input data matrix (generated by \code{\link{sim_create_x}})
#' @param covariates Integer vector that defines the types of covariates
#' (other than id and age). Different integers correspond to the
#' following covariate types:
#' \itemize{
#'   \item 0 = disease-related age
#'   \item 1 = other continuous covariate
#'   \item 2 = a categorical covariate that interacts with age
#'   \item 3 = a categorical covariate that acts as a group offset
#'   \item 4 = a categorical covariate that that acts as a group offset AND
#'   is restricted to have value 0 for controls and 1 for cases
#' }
#' @param relevances Relative relevance of each component. Must have be a vector
#' so that \cr
#'  \code{length(relevances) =  2 + length(covariates)}. \cr
#' First two values define the relevance of the individual-specific age and
#' shared age component, respectively.
#' @param lengthscales A vector so that \cr \code{length(lengthscales) = }
#' \code{2 + sum(covariates \%in\% c(0,1,2))}.
#' @param X_affected which individuals are affected by the disease
#' @param dis_fun A function or a string that defines the disease effect. If
#' this is a function, that function is used to generate the effect.
#' If \code{dis_fun} is "gp_vm" or "gp_ns", the disease component is drawn from
#' a nonstationary GP prior ("vm" is the variance masked version of it).
#' @param bin_kernel Should the binary kernel be used for categorical
#' covariates? If this is \code{TRUE}, the effect will exist only for group 1.
#' @param steepness Steepness of the input warping function. This is only used
#' if the disease component is in the model.
#' @param vm_params Parameters of the variance mask function. This is only
#' needed if \code{useMaskedVarianceKernel = TRUE}.
#' @param force_zeromean Should each component (excluding the disease age
#' component) be forced to have a zero mean?
#' @return a data frame FFF where one column corresponds to one additive
#' component
sim_create_f <- function(X,
                         covariates,
                         relevances,
                         lengthscales,
                         X_affected,
                         dis_fun,
                         bin_kernel,
                         steepness,
                         vm_params,
                         force_zeromean) {
  i4 <- which(covariates == 4)
  covariates[i4] <- 3
  D <- c(1, 2, covariates + 3)
  useMaskedVarianceKernel <- TRUE
  if (is.character(dis_fun)) {
    if (dis_fun == "gp_vm") {
      useMaskedVarianceKernel <- TRUE
    } else if (dis_fun == "gp_ns") {
      useMaskedVarianceKernel <- FALSE
    } else {
      stop("Invalid keyword for input dis_fun (", dis_fun, ")")
    }
  }

  KK <- sim_kernels(
    X, D, lengthscales,
    X_affected, bin_kernel,
    useMaskedVarianceKernel,
    steepness, vm_params
  )

  labs <- sim_name_components(D, colnames(X))
  FFF <- sim_draw_components(KK)
  if (sum(D == 3) == 1) {
    i_dis <- which(labs == "diseaseAge")
  } else {
    i_dis <- -1
  }
  if (i_dis > 0 && is.function(dis_fun)) {
    FFF[, i_dis] <- sim_disease_effect(X[, 1], X[, 3], dis_fun)
  } else {
    # do nothing, keep the component drawn from a GP
  }

  if (!bin_kernel) {
    i_zm <- c(which(D == 1), which(D == 5))
    i_skip <- c(i_dis, i_zm)
  } else {
    i_skip <- c(i_dis)
  }

  FFF <- sim_scale_relevances(FFF, relevances,
    force_zeromean = force_zeromean, i_skip
  )
  colnames(FFF) <- labs
  ret <- list(FFF = data.frame(FFF), KKK = KK)
  return(ret)
}


#' Compute all kernel matrices when simulating data
#'
#' @param X covariates
#' @param types vector of covariate types, so that
#' \itemize{
#'   \item 1 = ID
#'   \item 2 = age
#'   \item 3 = diseaseAge
#'   \item 4 = other continuous covariate
#'   \item 5 = a categorical covariate that interacts with age
#'   \item 6 = a categorical covariate that acts as an offset
#' }
#' @param lengthscales vector of lengthscales
#' @param X_affected which individuals are affected by the disease
#' @param bin_kernel whether or not binary (mask) kernel should be used for
#' categorical covariates
#' @param useMaskedVarianceKernel should the masked variance kernel be used
#' for drawing the disease component
#' @param steepness steepness of the input warping function
#' @param vm_params parameters of the variance mask function
#' @return a 3D array
sim_kernels <- function(X,
                        types,
                        lengthscales,
                        X_affected,
                        bin_kernel,
                        useMaskedVarianceKernel,
                        steepness,
                        vm_params) {
  n <- dim(X)[1]
  d <- length(types)
  KK <- array(0, c(n, n, d))
  t <- X[, which(types == 2)]
  id <- X[, which(types == 1)]
  n_ell <- sum(types != 6)
  d_ell <- length(lengthscales)
  if (d_ell != n_ell) {
    stop("lengthscales has length ", d_ell, ", should be ", n_ell)
  }
  j_ell <- 0
  ell <- lengthscales
  for (j in 1:d) {
    xj <- X[, j]
    if (types[j] == 1) {
      j_ell <- j_ell + 1
      N_tot <- length(unique(xj))
      Kj <- sim_kernel_zerosum(id, id, N_tot) *
        sim_kernel_se(t, t, ell = ell[j_ell])
    } else if (types[j] == 2) {
      j_ell <- j_ell + 1
      Kj <- sim_kernel_se(t, t, ell = ell[j_ell])
    } else if (types[j] == 3) {
      j_ell <- j_ell + 1
      Kj <- sim_kernel_bin(X_affected, X_affected) *
        sim_kernel_ns(xj, xj, ell = ell[j_ell], a = steepness, b = 0, c = 1)
      if (useMaskedVarianceKernel) {
        M <- sim_kernel_var_mask(xj, xj, vm_params, stp = steepness)
        Kj <- Kj * M
      }
    } else if (types[j] == 4) {
      j_ell <- j_ell + 1
      Kj <- sim_kernel_se(xj, xj, ell = ell[j_ell])
    } else if (types[j] == 5) {
      j_ell <- j_ell + 1
      if (bin_kernel) {
        Kj <- sim_kernel_bin(xj, xj) * sim_kernel_se(t, t, ell = ell[j_ell])
      } else {
        N_cat <- length(unique(xj))
        Kj <- sim_kernel_zerosum(xj, xj, N_cat) *
          sim_kernel_se(t, t, ell = ell[j_ell])
      }
    } else if (types[j] == 6) {
      if (bin_kernel) {
        Kj <- sim_kernel_bin(xj, xj)
      } else {
        N_cat <- length(unique(xj))
        Kj <- sim_kernel_zerosum(xj, xj, N_cat)
      }
    } else {
      stop("types contains invalid values")
    }
    KK[, , j] <- Kj
  }
  return(KK)
}

#' Create names for all components based on covariate names and types
#'
#' @param types vector of covariate types
#' @param names names of the covariates
#' @return a vector of component names
sim_name_components <- function(types, names) {
  d <- length(types)
  componentNames <- rep("foo", d)
  for (j in 1:d) {
    if (types[j] == 1) {
      lab <- "id*age"
    } else if (types[j] == 5) {
      lab <- paste(names[j], "*age", sep = "")
    } else {
      lab <- names[j]
    }
    componentNames[j] <- lab
  }
  return(componentNames)
}


#' Scale the effect sizes
#'
#' @param FFF matrix where one column corresponds to one additive component
#' @param relevances the desired variance of each component (column)
#' @param force_zeromean Should each component (excluding the disease age
#' component) be forced to have a zero mean.
#' @param i_skip induces of components for which the zero-mean forcing is
#' skipped
#' @return a new matrix \code{FFF}
sim_scale_relevances <- function(FFF, relevances,
                                 force_zeromean, i_skip) {

  # Some input checking
  d <- dim(FFF)[2]
  if (any(relevances < 0)) {
    stop("negative relevances not allowed!")
  }

  # Scale the columns to correct relevances
  for (j in 1:d) {
    std <- stats::sd(FFF[, j])
    if (std > 0) {
      FFF[, j] <- sqrt(relevances[j]) / std * FFF[, j]
      if (force_zeromean && !(j %in% i_skip)) {
        FFF[, j] <- FFF[, j] - mean(FFF[, j])
      }
    }
  }
  return(FFF)
}


#' Draw realizations of multivariate normals
#'
#' @param KK 3D matrix where \code{KK[,,j]} is the \code{j}th kernel matrix
#' @return a matrix \code{FFF}
sim_draw_components <- function(KK) {
  n <- dim(KK)[1]
  d <- dim(KK)[3]
  mu0 <- rep(0, n)
  FFF <- matrix(0, n, d)
  for (j in 1:d) {
    K <- KK[, , j]
    fj <- MASS::mvrnorm(1, mu0, K)
    FFF[, j] <- fj
  }
  return(FFF)
}

#' Draw disease component from a parameteric form
#'
#' @param X_id the id covariate
#' @param X_disAge the diseaseAge covariate
#' @param dis_fun the disease age effect function
#' @return a vector
sim_disease_effect <- function(X_id, X_disAge, dis_fun) {
  n <- length(X_id)
  uid <- unique(X_id)
  F_disAge <- rep(0, n)
  for (id in uid) {
    inds <- which(X_id == id)
    da <- X_disAge[inds]
    if (!is.nan(da[1])) {
      F_disAge[inds] <- dis_fun(da)
    }
  }
  return(F_disAge)
}

#' Simulate an input data frame X
#'
#' @param names Covariate names.
#' @inheritParams sim_create_x_D
#' @inheritParams sim_create_x_check
#' @inheritParams sim_create_x_dis_age
#' @inheritParams sim_create_x_other
#' @param t_jitter Standard deviation of the jitter added to the given
#' measurement times.
#' @return a list
sim_create_x <- function(N,
                         covariates,
                         names,
                         n_categs,
                         t_data,
                         t_jitter,
                         t_effect_range,
                         continuous_info) {

  # Validate input
  D <- sim_create_x_D(covariates)
  checked <- sim_create_x_check(n_categs, D, t_data, t_effect_range)

  # Initialize data
  k <- length(t_data)
  id <- rep(1:N, each = k)
  id <- as.factor(id)
  age <- sim_draw_measurement_times(N, t_data, t_jitter)
  X <- data.frame(id, age)

  # Effect times and disease ages
  parsed_dis <- sim_create_x_dis_age(X, k, N, D, age, checked$t_effect_range)
  dis_age <- parsed_dis$dis_age
  X <- parsed_dis$X

  # Other covariates
  cinfo <- continuous_info
  parsed_x <- sim_create_x_other(X, k, N, D, n_categs, dis_age, cinfo)
  X <- parsed_x$X
  colnames(X) <- names

  # Return
  list(
    info = checked$info,
    onsets = parsed_dis$teff,
    N_cases = parsed_dis$N_cases,
    par_cont = parsed_x$par_cont,
    X = X
  )
}

#' Helper function for sim_create_x
#'
#' @param covariates Integer vector that defines the types of covariates
#' (other than id and age). If not given, only the id and age
#' covariates are created. Different integers correspond to the following
#' covariate types:
#' \itemize{
#'   \item 0 = disease-related age
#'   \item 1 = other continuous covariate
#'   \item 2 = a categorical covariate that interacts with age
#'   \item 3 = a categorical covariate that acts as a group offset
#'   \item 4 = a categorical covariate that that acts as a group offset AND
#'   is restricted to have value 0 for controls and 1 for cases
#' }
#' @return an array of five integers
sim_create_x_D <- function(covariates) {
  D <- rep(0, 5)
  D[1] <- sum(covariates == 0)
  D[2] <- sum(covariates == 1)
  D[3] <- sum(covariates == 2)
  D[4] <- sum(covariates == 3)
  D[5] <- sum(covariates == 4)
  return(D)
}

#' Input check helper function for sim_create_x
#'
#' @param n_categs An integer vector defining the number of categories
#' for each categorical covariate, so that \code{length(n_categs)} equals to
#' the number of 2's and 3's in the \code{covariates} vector.
#' @param D covariate type array
#' @param t_data Measurement times.
#' @param t_effect_range Time interval from which the disease effect times are
#' sampled uniformly. Alternatively, This can any function that returns the
#' (possibly randomly generated) real disease effect time for one individual.
#' @return a list
sim_create_x_check <- function(n_categs, D, t_data, t_effect_range) {

  # Check length of n_categs
  L1 <- length(n_categs)
  L2 <- sum(D[3:4])
  if (L1 != L2) {
    msg <- paste0(
      "Length of <n_categs> must be same as the number of 2's",
      " and 3's in the <covariates> vector!",
      " Found = ", L1, ", should be = ", L2, "."
    )
    stop(msg)
  }

  # Check values of n_categs
  n_invalid <- sum(n_categs <= 1)
  if (n_invalid > 0) {
    stop("<n_categs> must only contain integers larger than 1!")
  }

  # Parse effect time range
  info <- ""
  if (is.character(t_effect_range)) {
    if (t_effect_range == "auto") {
      ran <- range(t_data)
      mrn <- mean(ran)
      t_effect_range <- c(mrn, mrn)
      if (D[1] == 1) {
        info <- paste0(
          "Disease effect time range set to [",
          t_effect_range[1], ", ", t_effect_range[2], "].\n"
        )
      }
    } else {
      stop("invalid t_effect_range!")
    }
  }

  # Return
  list(
    t_effect_range = t_effect_range,
    info = info
  )
}

#' Helper function for sim_create_x
#'
#' @param X current covariate matrix
#' @param k number of time points per individual
#' @param N Number of individuals.
#' @param D covariate type array
#' @param age the age covariate
#' @param t_effect_range Parsed version of \code{t_effect_range}
#' @return list
sim_create_x_dis_age <- function(X, k, N, D, age, t_effect_range) {
  if (D[1] > 0) {
    N_cases <- round(N / 2)
    if (is.function(t_effect_range)) {
      teff <- rep(0, N_cases)
      for (j in 1:N_cases) {
        teff[j] <- t_effect_range()
      }
    } else {
      teff <- stats::runif(N_cases,
        min = t_effect_range[1],
        max = t_effect_range[2]
      )
    }
    teff <- c(teff, rep(NaN, N - N_cases))
    dis_age <- sim_onsets_to_dis_age(teff, age, k)
    X <- cbind(X, dis_age)
  } else {
    dis_age <- NULL
    teff <- rep(NaN, N)
    N_cases <- NaN
  }

  # Return
  lt <- length(teff)
  names(teff) <- seq_len(lt)
  list(
    X = X,
    teff = teff,
    N_cases = N_cases,
    dis_age = dis_age
  )
}



#' Helper function for sim_create_x
#'
#' @inheritParams sim_create_x_dis_age
#' @inheritParams sim_create_x_check
#' @param continuous_info Info for generating continuous covariates. Must be a
#' list containing fields \code{lambda} and \code{mu}, which have length 3.
#' The continuous covariates are generated so that \code{x <- sin(a*t + b) + c},
#' where
#' \itemize{
#'   \item \code{t <- seq(0, 2*pi, length.out = k)}
#'   \item \code{a <- mu[1] + lambda[1]*stats::runif(1)}
#'   \item \code{b <- mu[2] + lambda[2]*stats::runif(1)}
#'   \item \code{c <- mu[3] + lambda[3]*stats::runif(1)}
#' }
#' @param dis_age a disease-age vector
#' @return a list
sim_create_x_other <- function(X, k, N, D, n_categs, dis_age, continuous_info) {
  if (D[2] > 0) {
    mu <- continuous_info$mu
    lambda <- continuous_info$lambda
    CONT <- sim_draw_continuous(N, k, D[2], mu, lambda)
    cont <- CONT$C
    par_cont <- list(a = CONT$A, b = CONT$B, offset = CONT$OFS)
    X <- cbind(X, cont)
  } else {
    par_cont <- list(a = NULL, b = NULL, offset = NULL)
  }

  if ((D[3] + D[4]) > 0) {
    categ <- sim_draw_categorical(N, k, n_categs)
    X <- cbind(X, categ)
  }

  if (D[5] > 0) {
    if (D[1] == 0) {
      stop("cannot include a 4 in covariates if 0 is not included!")
    }
    if (D[5] == 1) {
      group <- as.factor(as.numeric(!is.nan(dis_age)))
      X <- cbind(X, group)
    } else {
      stop("only one 4 can be included in the covariates vector")
    }
  }

  # Return
  list(X = X, par_cont = par_cont)
}


#'  Draw the age covariate
#'
#' @param N number of individuals
#' @param t_data a vector of length \code{k}
#' @param t_jitter Standard deviation of the jitter added to the given
#' measurement times.
#' @return a vector of length \code{N*k}
sim_draw_measurement_times <- function(N, t_data, t_jitter) {
  k <- length(t_data)
  age <- rep(0, N * k)
  if (t_jitter < 0) {
    stop("t_jitter must be positive!")
  }
  for (i in 1:N) {
    idx <- (1 + (i - 1) * k):(i * k)
    t <- t_data + stats::rnorm(k, mean = 0, sd = t_jitter)
    age[idx] <- sort(t)
  }
  return(age)
}


#' Compute the disease-related ages
#' @param age the age covariate, a vector of length \code{N*k}
#' @param onsets true disease effect times, a vector of length \code{N}
#' @param k number of measurements per individual
#' @return the diseaseAge covariate, a vector of length \code{N*k}
sim_onsets_to_dis_age <- function(onsets, age, k) {
  N <- length(onsets)
  n <- N * k
  diseaseAge <- rep(0, n)
  for (i in 1:N) {
    inds <- ((i - 1) * k + 1):(i * k)
    if (!is.nan(onsets[i])) {
      diseaseAge[inds] <- age[inds] - onsets[i]
    } else {
      diseaseAge[inds] <- NaN
    }
  }
  return(diseaseAge)
}


#' Indepedently draw continuous variables for each individual
#'
#' @param N number of individuals
#' @param k number of timepoints
#' @param D number of variables
#' @param mu a vector of length 3
#' @param lambda a vector of length 3
#' @return a matrix of size \code{N} x \code{D}
sim_draw_continuous <- function(N, k, D, mu, lambda) {
  C <- matrix(0, N * k, D)
  A <- matrix(0, N, D)
  B <- matrix(0, N, D)
  OFS <- matrix(0, N, D)
  for (j in 1:D) {
    for (i in c(1:N)) {
      inds <- ((i - 1) * k + 1):(i * k)
      ttt <- seq(0, 2 * pi, length.out = k)
      a <- mu[1] + lambda[1] * stats::runif(1)
      b <- mu[2] + lambda[2] * stats::runif(1)
      ofs <- mu[3] + lambda[3] * stats::runif(1)
      A[i, j] <- a
      B[i, j] <- b
      OFS[i, j] <- ofs
      C[inds, j] <- sin(a * ttt + b) + ofs
    }
  }
  return(list(A = A, B = B, OFS = OFS, C = C))
}


#' Independently draw categorical variables for each individual
#'
#' @param N number of individuals
#' @param k number of timepoints
#' @param v vector of numbers of different categories
#' @return a data frame with \code{N} rows and \code{D} = \code{length(v)}
#' columns
sim_draw_categorical <- function(N, k, v) {
  D <- length(v)
  C <- matrix(0, N, D)
  for (i in seq_len(D)) {
    C[, i] <- sample.int(v[i], N, replace = TRUE)
  }
  C_seq <- seq_len(nrow(C))
  rows <- rep(C_seq, each = k)
  C <- data.frame(C[rows, ])
  for (i in seq_len(D)) {
    C[, i] <- as.factor(C[, i])
  }
  return(C)
}

#' Compute a squared exponential kernel matrix
#'
#' @param x1 vector of length n
#' @param x2 vector of length m
#' @param alpha marginal std (default = 1)
#' @param ell lengthscale (default = 1)
#' @return A kernel matrix of size n x m
sim_kernel_se <- function(x1, x2, alpha = 1, ell = 1) {
  if (ell <= 0) {
    stop("ell must be positive!")
  }
  if (alpha < 0) {
    stop("alpha cannot be negative")
  }
  n1 <- length(x1)
  n2 <- length(x2)
  X1 <- matrix(rep(x1, each = n2), n1, n2, byrow = T)
  X2 <- matrix(rep(x2, n1), n1, n2, byrow = T)
  K <- alpha^2 * exp(-0.5 * (X1 - X2)^2 / ell^2)
  return(K)
}


#' Compute a zero-sum kernel matrix
#'
#' @param x1 (integer) vector of length n
#' @param x2 (integer) vector of length m
#' @param M number of categories
#' @param alpha marginal std (default = 1)
#' @return A (binary) kernel matrix of size n x m
sim_kernel_zerosum <- function(x1, x2, M, alpha = 1) {
  if (alpha < 0) {
    stop("alpha cannot be negative")
  }
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(0, n1, n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      if (x1[i] == x2[j]) {
        K[i, j] <- 1.0
      } else {
        K[i, j] <- -1.0 / (M - 1)
      }
    }
  }
  return(alpha^2 * K)
}

#' Compute a binary kernel matrix
#'
#' @param x1 (integer) vector of length n
#' @param x2 (integer) vector of length m
#' @param alpha marginal std (default = 1)
#' @param pos_class the positive class label
#' @return A kernel matrix of size n x m
sim_kernel_bin <- function(x1, x2 = NULL, alpha = 1, pos_class = 1) {
  if (alpha < 0) {
    stop("alpha cannot be negative")
  }
  n1 <- length(x1)
  n2 <- length(x2)
  X1 <- matrix(rep(x1, each = n2), n1, n2, byrow = T)
  X2 <- matrix(rep(x2, n1), n1, n2, byrow = T)
  K1 <- matrix(as.numeric(X1 == pos_class), n1, n2)
  K2 <- matrix(as.numeric(X2 == pos_class), n1, n2)
  return(alpha^2 * K1 * K2)
}

#' Compute a nonstationary kernel matrix using input warping
#'
#' @param x1 vector of length n
#' @param x2 vector of length m
#' @param alpha marginal std (default = 1)
#' @param ell lengthscale in the warped space
#' @param a steepness of the warping function rise
#' @param b location of the effective time window
#' @param c maximum range
#' @param nan_replace the value to use for replacing NaN values
#' @return A kernel matrix of size n x m
sim_kernel_ns <- function(x1, x2 = NULL, alpha = 1, ell, a, b, c,
                          nan_replace = 0) {
  if (a <= 0) {
    stop("a must be positive")
  }
  if (c <= 0) {
    stop("c must be positive")
  }
  x1[is.nan(x1)] <- nan_replace
  x2[is.nan(x2)] <- nan_replace
  w1 <- sim_warp_input(x1, a, b, c)
  w2 <- sim_warp_input(x2, a, b, c)
  K <- sim_kernel_se(w1, w2, alpha, ell)
  return(K)
}


#' Compute the multiplier matrix K_beta (to enable heterogeneous
#' disease effect)
#'
#' @param beta a row vector of length \code{N_cases}
#' @param row_to_caseID_1 mapping from row index to case ID
#' @param row_to_caseID_2 mapping from row index to case ID
#' @return a matrix
sim_kernel_beta <- function(beta, row_to_caseID_1, row_to_caseID_2) {
  n1 <- length(row_to_caseID_1)
  n2 <- length(row_to_caseID_2)
  BETA <- matrix(0, n1, n2)
  for (i in 1:n1) {
    i_case <- row_to_caseID_1[i]
    if (i_case > 0) {
      b1 <- beta[i_case]
    } else {
      b1 <- 0
    }
    for (j in 1:n2) {
      j_case <- row_to_caseID_2[j]
      if (j_case > 0) {
        b2 <- beta[j_case]
      } else {
        b2 <- 0
      }
      BETA[i, j] <- sqrt(b1 * b2)
    }
  }
  return(BETA)
}


#' Compute the variance mask kernel matrix
#'
#' @param disAge1 disease-related age covariate vector of length \code{n1}
#' @param disAge2 disease-related age covariate vector of length \code{n2}
#' @param vm_params vector of two mask function parameters
#' @param stp input warping steepness
#' @param nan_replace value to replace nans in disAge vectors
#' @return a matrix of size \code{n1} x \code{n2}
sim_kernel_var_mask <- function(disAge1, disAge2, vm_params,
                                stp, nan_replace = 0) {
  disAge1[is.nan(disAge1)] <- nan_replace
  disAge2[is.nan(disAge2)] <- nan_replace
  a <- stp * vm_params[2]
  h <- vm_params[1]
  if (h >= 1 || h <= 0) {
    stop("vm_params[1] must be between 0 and 1!")
  }
  if (vm_params[2] <= 0) {
    stop("vm_params[2] must be greater than 0!")
  }

  r <- 1 / a * log(h / (1 - h))
  s1 <- sim_var_mask(disAge1 - r, a)
  s2 <- sim_var_mask(disAge2 - r, a)
  M <- tcrossprod(s1, s2)
  return(M)
}


#' Input warping function
#'
#' @param t a vector
#' @param a steepness of the rise
#' @param b location of the effective time window
#' @param c maximum range
#' @return a vector of warped inputs \code{w(t)}
sim_warp_input <- function(t, a, b, c) {
  w <- 2 * c * (-0.5 + 1 / (1 + exp(-a * (t - b))))
  return(w)
}

#' Variance masking function
#'
#' @param x  vector of length \code{n}
#' @param a a positive real number
#' @return a vector of length \code{n}
sim_var_mask <- function(x, a) {
  y <- 1 / (1 + exp(-a * x))
  return(y)
}
