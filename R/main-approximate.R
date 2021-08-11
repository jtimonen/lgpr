#' Create an approximate GP model
#'
#' @export
#' @param model An object of class \linkS4class{lgpmodel}.
#' @param bf_options TODO!
#' @family main functions
#' @return An object of class \linkS4class{lgpmodel}.
approximate_model <- function(model, bf_options = list(M_bf = 30, L_bf = 3)) {
  si <- stan_input_approx(model, bf_options)
  model@stan_input <- si
  return(model)
}

# Create Stan data for approximate model
stan_input_approx <- function(model, bf_options) {
  si <- get_stan_input(model)
  comps <- dollar(si, "components")
  C_mats <- create_C_matrices(si)
  C_decs <- decompose_C_matrices(C_mats)
  si_add <- list(
    C_matrices = C_mats,
    C_eigvals = dollar(C_decs, "values"),
    C_eigvecs = dollar(C_decs, "vectors")
  )
  si <- c(si, si_add, bf_options)
  validate_approx_opts(bf_options)
  validate_stan_input_approx(si) # extra computation, could be commented
  return(si)
}

# Compute the C x C matrix for each term
create_C_matrices <- function(si) {
  STREAM <- get_stream()
  comp <- data.frame(dollar(si, "components"))
  x_cat <- dollar(si, "x_cat")
  x_cat_num_levels <- dollar(si, "x_cat_num_levels")
  C_matrices <- list()
  J <- nrow(comp)
  type <- dollar(comp, "type")
  ktype <- dollar(comp, "ker")
  icat <- dollar(comp, "cat")
  for (j in seq_len(J)) {
    idx <- icat[j]
    if (idx != 0) {
      z <- unique(x_cat[idx, ])
      M <- x_cat_num_levels[idx]
      C_j <- STAN_kernel_const(z, z, ktype[j], M, STREAM)
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
  for (j in seq_len(J)) {
    eg <- eigen(C_matrices[[j]])
    VALS[[j]] <- dollar(eg, "values")
    VECS[[j]] <- dollar(eg, "vectors")
  }
  list(values = VALS, vectors = VECS)
}

# Validation
validate_approx_opts <- function(bf_options) {
  L_bf <- dollar(bf_options, "L_bf")
  M_bf <- dollar(bf_options, "M_bf")
  check_positive(L_bf)
  check_positive(M_bf)
}

# Validation
validate_stan_input_approx <- function(si) {
  C_mats <- dollar(si, "C_matrices")
  C_vals <- dollar(si, "C_eigvals")
  C_vecs <- dollar(si, "C_eigvecs")
  J <- length(C_mats)
  for (j in seq_len(J)) {
    V <- C_vecs[[j]]
    D <- diag(C_vals[[j]])
    C_rec <- V %*% D %*% t(V)
    diff <- as.numeric(C_mats[[j]] - C_rec)
    mae <- max(abs(diff))
    if (mae > 1e-12) {
      msg <- paste0("Numerical problem in decomposition ", j, ", MAE=", mae)
      msg <- paste0(msg, ". Please report a bug! (validate_stan_input_approx)")
      stop(msg)
    }
  }
}
