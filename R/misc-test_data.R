#' Create data for unit tests of Stan kernel array functions
#'
#' @param N number of individuals
#' @return A list of data suitable for testing \code{STAN_kernel_all}
test_data_x <- function(N = 3) {
  x1_id <- rep(c(1, 2, 3), each = N)
  x1_z <- rep(c(0, 1, 1), each = N)
  x1_d <- rep(c(0, 1, 0), each = N)
  x1_disc <- list(x1_id, x1_z, x1_d)
  x1_cont <- list(rep(c(12, 24, 36), times = N))

  x2_id <- c(1, 1, 2, 3)
  x2_z <- c(0, 0, 1, 1)
  x2_d <- c(0, 0, 1, 0)
  x2_disc <- list(x2_id, x2_z, x2_d)
  x2_cont <- list(c(20, 20, 20, 20))

  num_levels <- c(N, 2, 2)
  covtypes <- c(0, 2, 2, 1, 3, 2)
  kertypes <- c(0, 0, 1, 0, 0, 2)
  cov_disc <- c(1, 2, 2, 0, 3, 2)
  cov_cont <- c(0, 1, 1, 1, 1, 1)
  comps <- list(covtypes, kertypes, cov_disc, cov_cont)

  dat <- list(
    x1_disc = x1_disc,
    x2_disc = x2_disc,
    x1_cont = x1_cont,
    x2_cont = x2_cont,
    num_levels = num_levels,
    components = comps
  )
  return(dat)
}

#' Create minimal data for Stan
#'
#' @return A list of data suitable for testing \code{lgp.stan}
minimal_stan_data <- function() {
  n <- 8
  array_dim0 <- function(x) {
    array(0, dim = c(0, x))
  }
  numeric_dim0 <- array(0, dim = c(0))

  stan_data <- list(
    is_verbose = 0,
    is_f_sampled = 0,
    is_heter = 0,
    is_uncrt = 0,
    is_likelihood_skipped = 0,
    is_vm_used = 0,
    is_generated_skipped = 0,
    num_obs = n,
    num_subjects = 1,
    num_cases = 0,
    num_cov_cont = 1,
    num_cov_disc = 0,
    num_comps = 1,
    num_ell = 1,
    num_dis = 0,
    obs_model = 1,
    components = matrix(c(1, 0, 0, 1), 4, 1, byrow = TRUE),
    y_cont = matrix(sin(0.5 * seq_len(n)), 1, n, byrow = TRUE),
    x_cont = matrix(seq_len(n), 1, n, byrow = TRUE),
    y_disc = array_dim0(n),
    x_disc = array_dim0(n),
    y_num_trials = array_dim0(n),
    num_levels = numeric_dim0,
    idx_expand = array_dim0(n),

    prior_alpha = matrix(c(1, 0), 1, 2, byrow = TRUE),
    prior_ell = matrix(c(6, 0), 1, 2, byrow = TRUE),
    prior_wrp = array_dim0(2),
    prior_sigma = matrix(c(6, 0), 1, 2, byrow = TRUE),
    prior_phi = array_dim0(2),
    prior_teff = array_dim0(3),

    hyper_alpha = matrix(c(0, 1, 0), 1, 3, byrow = TRUE),
    hyper_ell = matrix(c(1, 1, 0), 1, 3, byrow = TRUE),
    hyper_wrp = array_dim0(3),
    hyper_sigma = matrix(c(1, 1, 0), 1, 3, byrow = TRUE),
    hyper_phi = array_dim0(3),
    hyper_beta = array_dim0(2),
    hyper_teff = array_dim0(3),

    teff_obs = array_dim0(0),
    teff_lb = array_dim0(0),
    teff_ub = array_dim0(0),
    c_hat = rep(0, n),
    delta = 1e-8,
    vm_params = array_dim0(2)
  )
  return(stan_data)
}
