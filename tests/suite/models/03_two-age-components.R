# A model with two age components
library(lgpr)

# Create data using two age components
setup_data <- function(N = 10, H = 10, sigma = 0.5) {
  ttt <- seq(12, 72, by = H)
  K <- length(ttt)
  age <- rep(ttt, N) + rep(H / N * 1:N, each = K)
  id <- as.factor(rep(1:N, each = K))
  dat <- data.frame(id, age)
  f1 <- -0.001 * (-60 + age)**2 # slowly changing
  f2 <- exp(-sin(0.3 * age)) # quickly changing
  f1 <- f1 / sd(f1)
  f2 <- f2 / sd(f2)
  dat$age_fast <- age
  dat$y <- f1 + f2 + sigma * stats::rnorm(n = length(age))
  return(dat)
}

# Create model where there are different priors for the two lengthscales
setup_model <- function(dat, ...) {
  dat <- setup_data(...)
  d1 <- log_normal(1.0, 0.3)
  d2 <- normal(0, 0.3)
  prior <- list(ell = list(d1, d2))
  model <- create_model(y ~ age + age_fast, dat, prior = prior)
  return(model)
}

expected_relevances <- function() {
  c(0.0764848, 0.81305154, 0.1104637)
}
