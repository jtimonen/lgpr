# List available likelihood names
likelihood_list <- function() {
  c("gaussian", "poisson", "nb", "binomial", "bb")
}

# Convert Stan likelihood encoding to a string
likelihood_as_str <- function(index) {
  names <- likelihood_list()
  L <- length(names)
  check_interval(index, 1, L)
  name <- names[index]
  return(name)
}

# Convert likelihood name to Stan encoding
likelihood_as_int <- function(likelihood) {
  likelihood <- tolower(likelihood)
  allowed <- likelihood_list()
  index <- check_allowed(likelihood, allowed)
  return(index)
}

# Check if likelihood is binomial or beta-binomial
is_bin_or_bb <- function(likelihood) {
  likelihood %in% c("binomial", "bb")
}

# Check if likelihood is Poisson or negative binomial
is_pois_or_nb <- function(likelihood) {
  likelihood %in% c("poisson", "nb")
}

# Map x through link function
link <- function(x, likelihood) {
  check_allowed(likelihood, likelihood_list())
  if (is_pois_or_nb(likelihood)) {
    val <- log(x)
  } else if (is_bin_or_bb(likelihood)) {
    val <- log(x) - log(1 - x)
  } else {
    val <- x
  }
  return(val)
}

# Map x through inverse link function
link_inv <- function(x, likelihood) {
  allowed <- likelihood_list()
  check_allowed(likelihood, allowed)
  if (is_pois_or_nb(likelihood)) {
    x <- exp(x)
  } else if (is_bin_or_bb(likelihood)) {
    x <- 1 / (1 + exp(-x))
  }
  return(x)
}

# Map distribution of f to predictive distribution
map_f_to_y <- function(f_mean, f_std, sigma2, y_scl) {
  # Compute y_mean and y_std on normalized scale
  y_mean <- f_mean
  y_var <- add_to_columns(f_std^2, sigma2)
  y_std <- sqrt(y_var)

  # Scale y_mean and y_std to original scale and return
  list(
    mean = apply_scaling(y_scl, y_mean, inverse = TRUE),
    std = y_scl@scale * y_std
  )
}

# Map draws of f to prediction scale
map_f_to_h <- function(model, f, c_hat, reduce) {
  # Add c_hat
  num_draws <- dim(f)[1]
  f <- f + repvec(c_hat, num_draws)

  # Apply inverse link function and reduction
  likelihood <- get_obs_model(model)
  h <- link_inv(f, likelihood)
  h <- apply_reduce(h, reduce)
  return(h)
}

# Divide vector a elementwise by the vector of numbers of trials
divide_by_num_trials <- function(a, fit) {
  check_type(fit, "lgpfit")
  likelihood <- get_obs_model(fit)
  if (!is_bin_or_bb(likelihood)) {
    return(a)
  }
  y_num_trials <- get_num_trials(fit)
  check_lengths(a, y_num_trials)
  a / y_num_trials
}
