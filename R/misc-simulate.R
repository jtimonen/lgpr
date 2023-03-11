#' Generate an artificial longitudinal data set
#'
#' @export
#' @description Generate an artificial longitudinal data set.
#' @inheritParams sim.create_x
#' @inheritParams sim.create_f
#' @inheritParams sim.create_y
#' @param N_affected Number of diseased individuals that are affected by the
#' disease. This defaults to the number of diseased individuals. This argument
#' can only be given if \code{covariates} contains a zero.
#' @param t_observed Determines how the disease effect time is observed. This
#' can be any function that takes the real disease effect time as an argument
#' and returns the (possibly randomly generated) observed onset/initiation time.
#' Alternatively, this can be a string of the form \code{"after_n"} or
#' \code{"random_p"} or \code{"exact"}.
#' @param f_var variance of f
#' @param c_hat a constant added to f
#' @return An object of class \linkS4class{lgpsim}.
#' @examples
#' # Generate Gaussian data
#' dat <- simulate_data(N = 4, t_data = c(6, 12, 24, 36, 48), snr = 3)
#'
#' # Generate negative binomially (NB) distributed count data
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
                          gamma = 0.2,
                          N_affected = round(N / 2),
                          t_effect_range = "auto",
                          t_observed = "after_0",
                          c_hat = 0,
                          dis_fun = "gp_warp_vm",
                          bin_kernel = FALSE,
                          steepness = 0.5,
                          vm_params = c(0.025, 1),
                          continuous_info = list(
                            mu = c(pi / 8, pi, -0.5),
                            lambda = c(pi / 8, pi, 1)
                          ),
                          N_trials = 1,
                          force_zeromean = TRUE) {
  # Input checks
  noise_type <- tolower(noise_type)
  check_length_geq(t_data, 3)
  names <- sim.check_covariates(covariates, relevances, names, n_categs)
  bad <- N_affected > round(N / 2)
  if (bad) stop("N_affected cannot be greater than round(N/2)!")
  bad <- is.unsorted(covariates)
  if (bad) stop("The covariates vector must be increasing!")

  # Generate the covariates
  IN <- sim.create_x(
    N = N,
    covariates = covariates,
    names = names,
    n_categs = n_categs,
    t_data = t_data,
    t_jitter = t_jitter,
    t_effect_range = t_effect_range,
    continuous_info = continuous_info
  )

  # Compute X_affected
  k <- length(t_data)
  X_affected <- c(rep(1, N_affected * k), rep(0, (N - N_affected) * k))

  # Simulate the components F
  COMP <- sim.create_f(
    dollar(IN, "X"),
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
  FFF <- dollar(COMP, "FFF")

  # Total signal f is a (scaled) sum of the components plus a constant
  # - if there are zero components, define total signal as just noise
  f <- rowSums(FFF)
  SD <- stats::sd(f)
  f <- if (SD > 0) sqrt(f_var) / SD * f else stats::rnorm(n = length(f))
  f <- f + c_hat

  # Generate noise
  NOISY <- sim.create_y(noise_type, f,
    snr = snr, phi = phi, gamma = gamma,
    N_trials = N_trials
  )
  h <- dollar(NOISY, "h")
  y <- dollar(NOISY, "y")
  noise <- y - h

  # Create the output objects
  dat <- cbind(dollar(IN, "X"), y)
  rownames(dat) <- 1:(N * k)
  comp <- cbind(FFF, f, h, noise, y)

  # Real diseaseAges to observed ones
  OBSERVED <- sim.data_to_observed(dat, t_observed)

  # Signal and noise sums of squares
  SSR <- sum((h - mean(h))^2)
  SSE <- sum(noise^2)

  # Create rest of the fields
  teff_true <- dollar(IN, "onsets")
  teff_obs <- dollar(OBSERVED, "onsets_observed")
  teff <- list(true = teff_true, observed = teff_obs)

  info <- list(
    par_ell = lengthscales,
    par_cont = dollar(IN, "par_cont"),
    p_signal = SSR / (SSR + SSE),
    msg = dollar(IN, "info"),
    noise_type = noise_type
  )

  # Return S4 class object
  new("lgpsim",
    data = dollar(OBSERVED, "dat"),
    response = "y",
    components = comp,
    kernel_matrices = dollar(COMP, "KKK"),
    effect_times = teff,
    info = info
  )
}


#' Simulate noisy observations
#'
#' @param noise_type Either "gaussian", "poisson", "nb" (negative binomial),
#' "binomial", or "bb" (beta-binomial).
#' @param snr The desired signal-to-noise ratio. This argument is valid
#' only when \code{noise_type} is \code{"gaussian"}.
#' @param phi The inverse overdispersion parameter for negative binomial data.
#' The variance is \code{g + g^2/phi}.
#' @param gamma The dispersion parameter for beta-binomial data.
#' @param N_trials The number of trials parameter for binomial data.
#' @param f The underlying signal.
#' @return A list \code{out}, where
#' \itemize{
#'   \item \code{out$h} is \code{f} mapped through an inverse link function
#'   (times \code{N_trials} if \code{noise_type} is binomial or beta-binomial)
#'   \item \code{out$y} is the noisy response variable.
#' }
sim.create_y <- function(noise_type, f, snr, phi, gamma, N_trials) {
  L <- length(f)
  y <- rep(0, L)
  h <- link_inv(f, noise_type)
  if (noise_type == "gaussian") {
    sf <- stats::var(f)
    sigma_n <- 1 / sqrt(snr) * sqrt(sf)
    y_n <- stats::rnorm(n = L, mean = 0, sd = 1)
    y_n <- sigma_n * (y_n - mean(y_n)) / stats::sd(y_n)
    y <- h + y_n
  } else if (noise_type == "poisson") {
    y <- stats::rpois(n = L, lambda = h)
  } else if (noise_type == "nb") {
    y <- stats::rnbinom(n = L, mu = h, size = phi)
  } else if (noise_type == "binomial") {
    y <- stats::rbinom(n = L, size = N_trials, prob = h)
    h <- N_trials * h
  } else if (noise_type == "bb") {
    check_interval(gamma, 0, 1)
    gam_t <- (1 - gamma) / gamma
    alpha <- gam_t * h
    beta <- gam_t * (1.0 - h)
    prob <- stats::rbeta(n = L, alpha, beta)
    y <- stats::rbinom(n = L, size = N_trials, prob = prob)
    h <- N_trials * h
  } else {
    stop("Unknown noise_type! Please report a bug.")
  }
  NOISY <- list(h = h, y = y)
  return(NOISY)
}


#' Create an input data frame X for simulated data
#'
#' @param N Number of individuals.
#' @param t_data Measurement times (same for each individual, unless
#' \code{t_jitter > 0} in which case they are perturbed).
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
#' @param names Covariate names.
#' @param n_categs An integer vector defining the number of categories
#' for each categorical covariate, so that \code{length(n_categs)} equals to
#' the number of 2's and 3's in the \code{covariates} vector.
#' @param t_effect_range Time interval from which the disease effect times are
#' sampled uniformly. Alternatively, This can any function that returns the
#' (possibly randomly generated) real disease effect time for one individual.
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
#' @param t_jitter Standard deviation of the jitter added to the given
#' measurement times.
#' @return a list
sim.create_x <- function(N,
                         covariates,
                         names,
                         n_categs,
                         t_data,
                         t_jitter,
                         t_effect_range,
                         continuous_info) {
  # Validate input
  D <- sim.create_x_D(covariates)
  checked <- sim.create_x_check(n_categs, D, t_data, t_effect_range)

  # Initialize data
  k <- length(t_data)
  id <- rep(1:N, each = k)
  id <- as.factor(id)
  age <- sim.draw_measurement_times(N, t_data, t_jitter)
  X <- data.frame(id, age)

  # Effect times and disease ages
  et_range <- dollar(checked, "t_effect_range")
  parsed_dis <- sim.create_x_dis_age(X, k, N, D, age, et_range)
  dis_age <- parsed_dis$dis_age
  X <- parsed_dis$X

  # Other covariates
  cinfo <- continuous_info
  parsed_x <- sim.create_x_other(X, k, N, D, n_categs, dis_age, cinfo)
  X <- dollar(parsed_x, "X")
  colnames(X) <- names

  # Return
  list(
    info = dollar(checked, "info"),
    onsets = dollar(parsed_dis, "teff"),
    N_cases = dollar(parsed_dis, "N_cases"),
    par_cont = dollar(parsed_x, "par_cont"),
    X = X
  )
}

#' Simulate latent function components for longitudinal data analysis
#'
#' @param X input data matrix (generated by \code{\link{sim.create_x}})
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
sim.create_f <- function(X,
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
    if (dis_fun == "gp_warp_vm") {
      useMaskedVarianceKernel <- TRUE
    } else if (dis_fun == "gp_warp") {
      useMaskedVarianceKernel <- FALSE
    } else {
      msg <- paste0("dis_fun must be gp_warp or gp_warp_vm! found = ", dis_fun)
      stop(msg)
    }
  }

  KK <- sim.kernels(
    X, D, lengthscales,
    X_affected, bin_kernel,
    useMaskedVarianceKernel,
    steepness, vm_params
  )

  labs <- sim.name_components(D, colnames(X))
  FFF <- sim.draw_components(KK)
  if (sum(D == 3) == 1) {
    i_dis <- which(labs == "diseaseAge")
  } else {
    i_dis <- -1
  }
  if (i_dis > 0 && is.function(dis_fun)) {
    FFF[, i_dis] <- sim.disease_effect(X[, 1], X[, 3], dis_fun)
  } else {
    # do nothing, keep the component drawn from a GP
  }

  if (!bin_kernel) {
    i_zm <- c(which(D == 1), which(D == 5))
    i_skip <- c(i_dis, i_zm)
  } else {
    i_skip <- c(i_dis)
  }

  FFF <- sim.scale_relevances(FFF, relevances,
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
#' categorical covariates (if not, the zerosum kernel is used)
#' @param useMaskedVarianceKernel should the masked variance kernel be used
#' for drawing the disease component
#' @param steepness steepness of the input warping function
#' @param vm_params parameters of the variance mask function
#' @return a 3D array
sim.kernels <- function(X,
                        types,
                        lengthscales,
                        X_affected,
                        bin_kernel,
                        useMaskedVarianceKernel,
                        steepness,
                        vm_params) {
  pos_class <- 1 # bin kernel thing
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
      Kj <- kernel_zerosum(id, id, N_tot) * kernel_eq(t, t, ell = ell[j_ell])
    } else if (types[j] == 2) {
      j_ell <- j_ell + 1
      Kj <- kernel_eq(t, t, ell = ell[j_ell])
    } else if (types[j] == 3) {
      j_ell <- j_ell + 1
      ell_ns <- ell[j_ell]
      Kj <- kernel_bin(X_affected, X_affected, pos_class) *
        kernel_ns(xj, xj, ell = ell_ns, a = steepness)
      if (useMaskedVarianceKernel) {
        M <- kernel_varmask(xj, xj, steepness, vm_params)
        Kj <- Kj * M
      }
    } else if (types[j] == 4) {
      j_ell <- j_ell + 1
      Kj <- kernel_eq(xj, xj, ell = ell[j_ell])
    } else if (types[j] == 5) {
      j_ell <- j_ell + 1
      if (bin_kernel) {
        Kj <- kernel_bin(xj, xj, pos_class) * kernel_eq(t, t, ell = ell[j_ell])
      } else {
        N_cat <- length(unique(xj))
        Kj <- kernel_zerosum(xj, xj, N_cat) * kernel_eq(t, t, ell = ell[j_ell])
      }
    } else if (types[j] == 6) {
      if (bin_kernel) {
        Kj <- kernel_bin(xj, xj, pos_class)
      } else {
        N_cat <- length(unique(xj))
        Kj <- kernel_zerosum(xj, xj, N_cat)
      }
    } else {
      stop("types contains invalid values")
    }
    KK[, , j] <- Kj
  }
  return(KK)
}


# Input check for the covariate-related arguments of simulate_data
sim.check_covariates <- function(covariates, relevances, names, n_cat) {
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
    check_lengths(names, covariates)
    names <- c("id", "age", names)
  } else {
    # Generate default names
    names <- sim.generate_names(covariates)
  }
  check_non_negative_all(relevances)
  return(names)
}


# Generate names for covariates, given vector of covariate types
sim.generate_names <- function(covariates) {
  # Get default names
  names <- c("id", "age")
  def <- c("x", "z", "offset", "group")
  D <- sim.create_x_D(covariates)
  if (D[1] == 1) {
    names <- c(names, "diseaseAge")
  }
  if (D[2] > 0) {
    str <- paste0(def[1], 1:D[2])
    names <- if (D[2] > 1) c(names, str) else c(names, def[1])
  }
  if (D[3] > 0) {
    str <- paste0(def[2], 1:D[3])
    names <- if (D[3] > 1) c(names, str) else c(names, def[2])
  }
  if (D[4] > 0) {
    str <- paste0(def[3], 1:D[4])
    names <- if (D[4] > 1) c(names, str) else c(names, def[3])
  }
  if (D[5] > 0) {
    names <- c(names, "group")
  }

  return(names)
}

# Convert real generated disease ages to observed ones
sim.data_to_observed <- function(dat, t_observed) {
  flag <- !("diseaseAge" %in% colnames(dat))
  id <- dollar(dat, "id")
  uid <- unique(id)
  N <- length(uid)
  onsets_observed <- rep(NaN, N)
  names(onsets_observed) <- c(1:N)
  if (flag) {
    # No modifications to data frame needed
    ret <- list(dat = dat, onsets_observed = onsets_observed)
    return(ret)
  } else {
    age <- dollar(dat, "age")
    dis_age <- dollar(dat, "diseaseAge")
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
          oo <- sim.apply_obs_onset_fun(t_observed, t0_real, rem)
          t0 <- dollar(oo, "t0")
          idx0 <- dollar(oo, "idx0")
        } else {
          parsed <- sim.parse_t_obs(t_observed)
          typ <- dollar(parsed, "type")
          val <- dollar(parsed, "value")
          if (typ == "random") {
            idx0 <- sim.rtgeom(length(irem), val)
            sim.check_too_far(idx0, rem)
            t0 <- rem[idx0]
          } else if (typ == "after") {
            idx0 <- val + 1
            sim.check_too_far(idx0, rem)
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

# Check if trying to simulate onset happening after measurement interval
#
# @param fun_obs Function for generating the observed onset
# from real onset.
# @param t0_real The real onset.
# @param rem Measurement points after real onset.
# @return A list with names \code{idx0} and \code{t0}.
sim.apply_obs_onset_fun <- function(fun_obs, t0_real, rem) {
  t_possible <- fun_obs(t0_real)
  inds0 <- which(rem > t_possible)
  if (length(inds0) < 1) {
    stop("There are no data points after t = ", t_possible, "!")
  } else {
    idx0 <- inds0[1] # next meas. point after detection is possible
    t0 <- rem[idx0]
  }
  list(idx0 = idx0, t0 = t0)
}

# Check if trying to simulate onset happening after measurement interval
#
# @param idx index of onset point relative to real onset
# @param rem measurement points after real onset
sim.check_too_far <- function(idx, rem) {
  check_length(idx, 1)
  if (idx > length(rem)) {
    stop("Not enough data points to go that far!")
  }
  TRUE
}

# Sample from the 'truncated geometric' distribution
# @param p a number between 0 and 1
# @param s an integer
# @param n number of samples
# @return an integer from the interval 1...n
sim.rtgeom <- function(s, p, n = 1) {
  s_seq <- 0:(s - 1)
  prob <- p^s_seq
  r <- sample.int(n = s, size = n, prob = prob, replace = TRUE)
  return(r)
}


# Parse the t_observed argument of simulate_data
sim.parse_t_obs <- function(t_observed) {
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



# Create names for all components based on covariate names and types
sim.name_components <- function(types, names) {
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


# Scale the effect sizes
#
# @param FFF Matrix where one column corresponds to one additive component
# @param relevances The desired variance of each component (column)
# @param force_zeromean Should each component (excluding the disease age
# component) be forced to have a zero mean?
# @param i_skip Indices of components for which the zero-mean forcing is
# skipped
# @return a new matrix like FFF
sim.scale_relevances <- function(FFF, relevances,
                                 force_zeromean, i_skip) {
  # Some input checking
  d <- dim(FFF)[2]
  check_non_negative_all(relevances)

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


# Draw realizations of multivariate normals
#
# @param KK 3D matrix where \code{KK[,,j]} is the \code{j}th kernel matrix
# @return A matrix FFF
sim.draw_components <- function(KK) {
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

# Draw disease component from a parameteric form (dis_fun)
sim.disease_effect <- function(X_id, X_disAge, dis_fun) {
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

# Helper function for sim.create_x
sim.create_x_D <- function(covariates) {
  D <- rep(0, 5)
  D[1] <- sum(covariates == 0)
  D[2] <- sum(covariates == 1)
  D[3] <- sum(covariates == 2)
  D[4] <- sum(covariates == 3)
  D[5] <- sum(covariates == 4)
  return(D)
}

# Input check helper function for sim.create_x (D is covariate type array)
sim.create_x_check <- function(n_categs, D, t_data, t_effect_range) {
  # Check length of n_categs
  L1 <- length(n_categs)
  L2 <- sum(D[3:4])
  msg <- paste0(
    "Length of <n_categs> must be same as the number of 2's",
    " and 3's in the <covariates> vector!",
    " Found = ", L1, ", should be = ", L2, "."
  )
  if (L1 != L2) {
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

# Helper function for sim.create_x
#
# @param X current covariate matrix
# @param k number of time points per individual
# @param N Number of individuals.
# @param D covariate type array
# @param age the age covariate
# @param t_effect_range Parsed version of \code{t_effect_range}
# @return list
sim.create_x_dis_age <- function(X, k, N, D, age, t_effect_range) {
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
    dis_age <- sim.onsets_to_dis_age(teff, age, k)
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



# Helper function for sim.create_x
#
# @inheritParams sim.create_x_dis_age
# @inheritParams sim.create_x_check
# @param dis_age a disease-age vector
# @return a list
sim.create_x_other <- function(X, k, N, D, n_categs, dis_age, continuous_info) {
  if (D[2] > 0) {
    mu <- dollar(continuous_info, "mu")
    lambda <- dollar(continuous_info, "lambda")
    CONT <- sim.draw_continuous(N, k, D[2], mu, lambda)
    cont <- dollar(CONT, "C")
    par_cont <- list(
      a = dollar(CONT, "A"),
      b = dollar(CONT, "B"),
      offset = dollar(CONT, "OFS")
    )
    X <- cbind(X, cont)
  } else {
    par_cont <- list(a = NULL, b = NULL, offset = NULL)
  }

  if ((D[3] + D[4]) > 0) {
    categ <- sim.draw_categorical(N, k, n_categs)
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


# Draw the age covariate
#
# @param N Number of individuals
# @param t_data A vector of length \code{k}
# @param t_jitter Standard deviation of the jitter added to the given
# measurement times.
# @return A vector of length \code{N*k}
sim.draw_measurement_times <- function(N, t_data, t_jitter) {
  k <- length(t_data)
  age <- rep(0, N * k)
  check_non_negative(t_jitter)
  for (i in 1:N) {
    idx <- (1 + (i - 1) * k):(i * k)
    t <- t_data + stats::rnorm(k, mean = 0, sd = t_jitter)
    age[idx] <- sort(t)
  }
  return(age)
}


# Compute the disease-related ages
# @param age The age covariate, a vector of length \code{N*k}
# @param onsets True disease effect times, a vector of length \code{N}
# @param k Number of measurements per individual
# @return The diseaseAge covariate, a vector of length \code{N*k}
sim.onsets_to_dis_age <- function(onsets, age, k) {
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


# Indepedently draw continuous variables for each individual
#
# @param N Number of individuals
# @param k Number of timepoints
# @param D Number of variables
# @param mu A vector of length 3
# @param lambda A vector of length 3
# @return A matrix of size \code{N} x \code{D}
sim.draw_continuous <- function(N, k, D, mu, lambda) {
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


# Independently draw categorical variables for each individual
#
# @param N Number of individuals
# @param k Number of timepoints
# @param v Vector of numbers of different categories
# @return A data frame with \code{N} rows and \code{D} = \code{length(v)}
# columns
sim.draw_categorical <- function(N, k, v) {
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
