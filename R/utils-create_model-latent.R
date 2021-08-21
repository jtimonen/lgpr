# Create a latent GP model from base model
create_model.latent <- function(base_model, likelihood, prior,
                                c_hat, num_trials, approx, verbose) {
  stan_input <- list(approx = approx)
  si <- get_stan_input(base_model)

  # Base latent model
  parsed <- standata_latent(base_model, likelihood, prior, c_hat, num_trials)
  si_add <- dollar(parsed, "to_stan")
  si <- c(si, si_add)
  model <- new("LatentGPModel",
    base_model,
    parsed_input = si,
    y_scaling = dollar(parsed, "y_scaling"),
    info = creation_info()
  )

  # Approximate model
  if (!is.null(approx)) {
    model <- create_model.latent_approx(model, approx)
  }
  return(model)
}

# Parse Stan data for latent model
standata_latent <- function(base_model, likelihood, prior, c_hat, num_trials) {
  dat <- get_data(base_model)
  y_name <- get_y_name(base_model)
  y_info <- standata_latent.likelihood(
    dat, likelihood, c_hat, num_trials, y_name
  )
  return(y_info)
}

# Parse the response variable and its likelihood model
standata_latent.likelihood <- function(data, likelihood, c_hat, num_trials,
                                       y_name) {
  LH <- likelihood_as_int(likelihood)

  # Check that data contains the response variable,
  # which is numeric and compatible with observation model
  check_in_data(y_name, data, "data")
  Y_RAW <- dollar(data, y_name)
  check_response(Y_RAW, LH)

  # Parse y and related inputs
  parse_y_latent(Y_RAW, y_name, LH, c_hat, num_trials)
}


# Parse raw response taken from input data frame (latent GP model)
parse_y_latent <- function(Y_RAW, y_name, LH, c_hat, num_trials) {
  N <- length(Y_RAW)
  if (LH == 1) {
    y_int <- array(Y_RAW, dim = c(0, N))
    y_real <- array(Y_RAW, dim = c(1, N))
  } else {
    y_int <- array(Y_RAW, dim = c(1, N))
    y_real <- array(Y_RAW, dim = c(0, N))
  }
  num_trials <- set_num_trials(num_trials, Y_RAW, LH)
  c_hat <- set_c_hat(c_hat, Y_RAW, LH, num_trials)

  # Create stan input parts
  to_stan <- list(
    N = N,
    y_int = y_int,
    y_real = y_real,
    obs_model = LH,
    y_num_trials = num_trials,
    c_hat = c_hat
  )

  # Return stan input and dummy normalizer, which is identity mapping
  list(to_stan = to_stan, y_scaling = new("lgpscaling", var_name = y_name))
}

# Convert given c_hat input to Stan input format
set_c_hat <- function(c_hat, response, LH, num_trials) {
  nb_or_pois <- LH %in% c(2, 3)
  binomial <- LH %in% c(4, 5)
  gaussian <- LH == 1

  # Create  c_hat if it is NULL
  if (is.null(c_hat)) {
    if (nb_or_pois) {
      c_hat <- log(mean(response)) # Poisson or  NB
    } else if (binomial) {
      check_all_leq(response, num_trials)
      p <- mean(response / num_trials)
      c_hat <- log(p / (1 - p)) # Binomial or BB
    } else {
      stopifnot(gaussian)
      c_hat <- mean(response) # Gaussian
    }
  }
  n <- length(response)
  L <- length(c_hat)

  # Throw error if given c_hat is misspecified
  if (L != 1 && L != n) {
    stop(
      "Invalid length of <c_hat>! Must be 1 or equal to number ",
      "of observartions (", n, "). Found = ", L
    )
  }
  if (L == 1) c_hat <- rep(c_hat, n) # given c_hat is one number
  return(c_hat)
}

# Convert given num_trials input to Stan input format
set_num_trials <- function(num_trials, y, LH) {
  N <- length(y)
  is_binom <- LH %in% c(4, 5)
  if (is.null(num_trials)) {
    num_trials <- rep(1, N)
  } else {
    if (!is_binom) {
      msg <- paste0(
        "Only give the <num_trials> argument if likelihood is",
        " 'binomial' or 'bb' (beta-binomial)!"
      )
      stop(msg)
    }
    L <- length(num_trials)
    if (L == 1) {
      num_trials <- rep(num_trials, N)
    } else if (L != N) {
      stop(
        "Invalid length of <num_trials>! Must be 1 or equal to number ",
        "of observartions (", N, "). Found = ", L
      )
    }
  }
  DIM <- as.numeric(is_binom)
  num_trials <- array(num_trials, dim = c(DIM, N))
  return(num_trials)
}


# Create a latent approximate GP model from latent GP model
create_model.latent_approx <- function(latent_model, approx) {
  si <- get_stan_input(latent_model)
  si_add <- approx # TODO: parse
  si <- c(si, si_add)
  new("LatentGPModelApprox",
    latent_model,
    parsed_input = si,
    info = creation_info()
  )
}
