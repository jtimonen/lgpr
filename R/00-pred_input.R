#' Create input for Stan model that computes predictions
#'
#' @inheritParams pred
#' @return a list that can be used as data for \code{lgp_predict.stan}
#' @name pred_input
NULL

#' @rdname pred_input
pred_input <- function(fit, x, reduce, draws) {
  si <- get_stan_input(fit)
  si_x_pred <- pred_input_x(fit, x)
  si_draws <- pred_input_draws(fit, reduce, draws)
  c(si, si_x_pred, si_draws)
}

#' @rdname pred_input
pred_input_x <- function(fit, x) {
  si <- get_stan_input(fit)
  if (is.null(x)) {
    out <- list(
      num_pred = dollar(si, "num_obs"),
      x_cont_PRED = dollar(si, "x_cont"),
      x_cont_unnorm_PRED = dollar(si, "x_cont_unnorm"),
      x_cont_mask_PRED = dollar(si, "x_cont_mask"),
      x_cat_PRED = dollar(si, "x_cat"),
      idx_expand_PRED = dollar(si, "idx_expand")
    )
  } else {
    num_pred <- nrow(x)
    stop("create_pred_stan_input.x_pred not implemented when x != NULL")
  }
  return(out)
}

#' @rdname pred_input
pred_input_draws <- function(fit, reduce, draws) {

  # Get param sets
  d_common <- pred_input_draws.common(fit, reduce, draws)
  if (is_f_sampled(fit)) {
    d_add <- pred_input_draws.latent(fit, reduce, draws)
  } else {
    d_add <- pred_input_draws.marginal(fit, reduce, draws)
  }
  c(d_common, d_add)
}

#' @rdname pred_input
pred_input_draws.marginal <- function(fit, reduce, draws) {
  S <- get_num_paramsets(fit, draws, reduce)
  list(
    d_sigma = get_draw_arr(fit, draws, reduce, "sigma", S, 1)
  )
}

#' @rdname pred_input
pred_input_draws.latent <- function(fit, reduce, draws) {
  S <- get_num_paramsets(fit, draws, reduce)
  si <- get_stan_input(fit)
  LH <- dollar(si, "obs_model")
  num_sigma <- as.numeric(LH == 1)
  num_phi <- as.numeric(LH == 3)
  num_gamma <- as.numeric(LH == 5)
  num_f <- 1


  # Get draws
  return(NULL)
}

#' @rdname pred_input
pred_input_draws.common <- function(fit, reduce, draws) {

  # Get dimensions
  S <- get_num_paramsets(fit, draws, reduce)
  si <- get_stan_input(fit)
  num_comps <- dollar(si, "num_comps")
  num_ell <- dollar(si, "num_ell")
  num_ns <- dollar(si, "num_ns")
  UNCRT <- dollar(si, "num_uncrt") > 0
  HETER <- dollar(si, "num_heter") > 0
  num_bt <- dollar(si, "num_bt")

  # Get draws
  list(
    S = S,
    d_alpha = get_draw_arr(fit, draws, reduce, "alpha", S, num_comps),
    d_ell = get_draw_arr(fit, draws, reduce, "ell", S, num_ell),
    d_wrp = get_draw_arr(fit, draws, reduce, "wrp", S, num_ns),
    d_beta = get_draw_arr_vec(fit, draws, reduce, "beta", S, HETER, num_bt),
    d_teff = get_draw_arr_vec(fit, draws, reduce, "teff", S, UNCRT, num_bt)
  )
}

#' @rdname pred_input
get_draw_arr <- function(fit, draws, reduce, par_name, S, D) {
  out <- array(0.0, dim = c(S, D))
  if (D > 0) {
    out <- get_draws(fit, draws = draws, reduce = reduce, pars = par_name)
  }
  return(out)
}

#' @rdname pred_input
get_draw_arr_vec <- function(fit, draws, reduce, par_name, S, B, V) {
  D <- as.numeric(B)
  out <- array(0.0, dim = c(S, D, V))
  if (B) {
    tmp <- get_draws(fit, draws = draws, reduce = reduce, pars = par_name)
    # tmp has dim c(S, V)
    out[, 1, ] <- tmp
  }
  return(out)
}
