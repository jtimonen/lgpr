
# Create additional Stan input for approximate models
stan_input_approx_precomp <- function(stan_input) {
  num_bf <- dollar(stan_input, "num_bf")
  if (any(num_bf > 0)) {
    comps <- dollar(stan_input, "components")

    # Categorical decomposition stuff
    C_mats <- create_C_matrices(stan_input)
    C_decs <- decompose_C_matrices(C_mats)
    validate_stan_input_approx(C_decs, C_mats) # extra computation
    si_add <- vectorize_C_decs(C_decs)

    # Basis function stuff
    J <- nrow(comps)
    num_xi <- num_bf * J
    si_add <- c(si_add, list(num_xi = num_xi))
  } else {
    si_add <- c()
  }
  return(si_add)
}

# Compute the C x C matrix for each term
create_C_matrices <- function(si) {
  STREAM <- get_stream()
  comp <- dollar(si, "components")
  Z <- dollar(si, "Z")
  Z_M <- dollar(si, "Z_M")
  C_matrices <- list()
  J <- nrow(comp)
  iz <- comp[, 1]
  kz <- comp[, 2]
  for (j in seq_len(J)) {
    idx <- iz[j]
    if (idx != 0) {
      z <- unique(Z[idx, ])
      M <- Z_M[idx]
      C_j <- STAN_kernel_const(z, z, kz[j], M, STREAM)
    } else {
      C_j <- matrix(1.0, 1, 1)
    }
    C_matrices[[j]] <- C_j
  }
  return(C_matrices)
}

# Compute the eigendecompositions for the C x C matrices for each term
decompose_C_matrices <- function(C_matrices) {
  J <- length(C_matrices)
  VALS <- list()
  VECS <- list()
  SIZES <- rep(1, J)
  for (j in seq_len(J)) {
    eg <- eigen(C_matrices[[j]])
    VALS[[j]] <- dollar(eg, "values")
    VECS[[j]] <- dollar(eg, "vectors")
    SIZES[j] <- length(VALS[[j]])
  }
  list(values = VALS, vectors = VECS, sizes = SIZES)
}


# Vectorize
vectorize_C_decs <- function(C_decs) {
  D <- dollar(C_decs, "values")
  V <- dollar(C_decs, "vectors")
  C <- dollar(C_decs, "sizes")
  J <- length(C)
  C_eigvals <- c()
  C_eigvecs <- c()
  for (j in seq_len(J)) {
    V_j <- as.vector(V[[j]]) # takes the matrix by columns
    C_eigvals <- c(C_eigvals, D[[j]])
    C_eigvecs <- c(C_eigvecs, V_j)
  }
  list(
    C_eigvals = C_eigvals,
    C_eigvecs = C_eigvecs,
    C_sizes = dollar(C_decs, "sizes"),
    len_eigvals = length(C_eigvals),
    len_eigvecs = length(C_eigvecs)
  )
}


# Validation
validate_stan_input_approx <- function(C_decs, C_mats) {
  C_vals <- dollar(C_decs, "values")
  C_vecs <- dollar(C_decs, "vectors")
  J <- length(C_mats)
  for (j in seq_len(J)) {
    V <- C_vecs[[j]]
    D <- diag(C_vals[[j]])
    C_rec <- V %*% D %*% t(V)
    diff <- as.vector(C_mats[[j]] - C_rec)
    mae <- max(abs(diff))
    if (mae > 1e-12) {
      msg <- paste0("Numerical problem in decomposition ", j, ", MAE=", mae)
      msg <- paste0(msg, ". Please report a bug! (validate_stan_input_approx)")
      stop(msg)
    }
  }
}
