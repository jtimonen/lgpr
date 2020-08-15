// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "lgpr_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// STAN_check_real_positive
void STAN_check_real_positive(const double& a, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_check_real_positive(SEXP aSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    STAN_check_real_positive(a, pstream__);
    return R_NilValue;
END_RCPP
}
// STAN_check_prob_positive
void STAN_check_prob_positive(const double& a, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_check_prob_positive(SEXP aSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    STAN_check_prob_positive(a, pstream__);
    return R_NilValue;
END_RCPP
}
// STAN_warp_input
Eigen::Matrix<double, Eigen::Dynamic, 1> STAN_warp_input(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x, const double& a, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_warp_input(SEXP xSEXP, SEXP aSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_warp_input(x, a, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_var_mask
Eigen::Matrix<double, Eigen::Dynamic, 1> STAN_var_mask(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x, const double& a, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_var_mask(SEXP xSEXP, SEXP aSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_var_mask(x, a, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_expand
Eigen::Matrix<double, Eigen::Dynamic, 1> STAN_expand(const Eigen::Matrix<double, Eigen::Dynamic, 1>& v, const std::vector<int>& idx_expand, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_expand(SEXP vSEXP, SEXP idx_expandSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type idx_expand(idx_expandSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_expand(v, idx_expand, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_edit_x_cont
Eigen::Matrix<double, Eigen::Dynamic, 1> STAN_edit_x_cont(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x_cont, const std::vector<int>& idx_expand, const Eigen::Matrix<double, Eigen::Dynamic, 1>& teff_obs, const Eigen::Matrix<double, Eigen::Dynamic, 1>& teff, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_edit_x_cont(SEXP x_contSEXP, SEXP idx_expandSEXP, SEXP teff_obsSEXP, SEXP teffSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type x_cont(x_contSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type idx_expand(idx_expandSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type teff_obs(teff_obsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type teff(teffSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_edit_x_cont(x_cont, idx_expand, teff_obs, teff, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_kernel_base_zerosum
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> STAN_kernel_base_zerosum(const std::vector<int>& x1, const std::vector<int>& x2, const int& num_cat, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_kernel_base_zerosum(SEXP x1SEXP, SEXP x2SEXP, SEXP num_catSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const int& >::type num_cat(num_catSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_kernel_base_zerosum(x1, x2, num_cat, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_kernel_base_cat
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> STAN_kernel_base_cat(const std::vector<int>& x1, const std::vector<int>& x2, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_kernel_base_cat(SEXP x1SEXP, SEXP x2SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_kernel_base_cat(x1, x2, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_kernel_base_bin_mask
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> STAN_kernel_base_bin_mask(const std::vector<int>& x1, const std::vector<int>& x2, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_kernel_base_bin_mask(SEXP x1SEXP, SEXP x2SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_kernel_base_bin_mask(x1, x2, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_kernel_const
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> STAN_kernel_const(const std::vector<int>& x1, const std::vector<int>& x2, const int& kernel_type, const int& num_cat, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_kernel_const(SEXP x1SEXP, SEXP x2SEXP, SEXP kernel_typeSEXP, SEXP num_catSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const int& >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_cat(num_catSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_kernel_const(x1, x2, kernel_type, num_cat, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_kernel_const_all
std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > STAN_kernel_const_all(const int& n1, const int& n2, const std::vector<std::vector<int> >& x1, const std::vector<std::vector<int> >& x2, const std::vector<std::vector<int> >& x1_mask, const std::vector<std::vector<int> >& x2_mask, const std::vector<int>& num_levels, const std::vector<std::vector<int> >& components, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_kernel_const_all(SEXP n1SEXP, SEXP n2SEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP x1_maskSEXP, SEXP x2_maskSEXP, SEXP num_levelsSEXP, SEXP componentsSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const int& >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int> >& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int> >& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int> >& >::type x1_mask(x1_maskSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int> >& >::type x2_mask(x2_maskSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type num_levels(num_levelsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int> >& >::type components(componentsSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_kernel_const_all(n1, n2, x1, x2, x1_mask, x2_mask, num_levels, components, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_kernel_base_var_mask
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> STAN_kernel_base_var_mask(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x1, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x2, const double& steepness, const std::vector<double>& vm_params, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_kernel_base_var_mask(SEXP x1SEXP, SEXP x2SEXP, SEXP steepnessSEXP, SEXP vm_paramsSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const double& >::type steepness(steepnessSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type vm_params(vm_paramsSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_kernel_base_var_mask(x1, x2, steepness, vm_params, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_kernel_all
std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > STAN_kernel_all(const int& n1, const int& n2, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& K_const, const std::vector<std::vector<int> >& components, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& x1, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& x2, const std::vector<double>& alpha, const std::vector<double>& ell, const std::vector<double>& wrp, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& beta, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& teff, const std::vector<std::vector<double> >& vm_params, const std::vector<int>& idx1_expand, const std::vector<int>& idx2_expand, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& teff_obs, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_kernel_all(SEXP n1SEXP, SEXP n2SEXP, SEXP K_constSEXP, SEXP componentsSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP alphaSEXP, SEXP ellSEXP, SEXP wrpSEXP, SEXP betaSEXP, SEXP teffSEXP, SEXP vm_paramsSEXP, SEXP idx1_expandSEXP, SEXP idx2_expandSEXP, SEXP teff_obsSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const int& >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& >::type K_const(K_constSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int> >& >::type components(componentsSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type wrp(wrpSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& >::type teff(teffSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double> >& >::type vm_params(vm_paramsSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type idx1_expand(idx1_expandSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type idx2_expand(idx2_expandSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& >::type teff_obs(teff_obsSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_kernel_all(n1, n2, K_const, components, x1, x2, alpha, ell, wrp, beta, teff, vm_params, idx1_expand, idx2_expand, teff_obs, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_gp_posterior_helper
std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > STAN_gp_posterior_helper(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& Ly, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& K, const Eigen::Matrix<double, Eigen::Dynamic, 1>& v, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_gp_posterior_helper(SEXP LySEXP, SEXP KSEXP, SEXP vSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type Ly(LySEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_gp_posterior_helper(Ly, K, v, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_gp_posterior
std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > STAN_gp_posterior(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& KX, const Eigen::Matrix<double, Eigen::Dynamic, 1>& y, const double& delta, const double& sigma, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_gp_posterior(SEXP KXSEXP, SEXP ySEXP, SEXP deltaSEXP, SEXP sigmaSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& >::type KX(KXSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_gp_posterior(KX, y, delta, sigma, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// STAN_log_prior
double STAN_log_prior(const double& x, const std::vector<int>& types, const std::vector<double>& hyper, std::ostream* pstream__);
RcppExport SEXP _lgpr_STAN_log_prior(SEXP xSEXP, SEXP typesSEXP, SEXP hyperSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type types(typesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type hyper(hyperSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(STAN_log_prior(x, types, hyper, pstream__));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4lgp_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_lgpr_STAN_check_real_positive", (DL_FUNC) &_lgpr_STAN_check_real_positive, 2},
    {"_lgpr_STAN_check_prob_positive", (DL_FUNC) &_lgpr_STAN_check_prob_positive, 2},
    {"_lgpr_STAN_warp_input", (DL_FUNC) &_lgpr_STAN_warp_input, 3},
    {"_lgpr_STAN_var_mask", (DL_FUNC) &_lgpr_STAN_var_mask, 3},
    {"_lgpr_STAN_expand", (DL_FUNC) &_lgpr_STAN_expand, 3},
    {"_lgpr_STAN_edit_x_cont", (DL_FUNC) &_lgpr_STAN_edit_x_cont, 5},
    {"_lgpr_STAN_kernel_base_zerosum", (DL_FUNC) &_lgpr_STAN_kernel_base_zerosum, 4},
    {"_lgpr_STAN_kernel_base_cat", (DL_FUNC) &_lgpr_STAN_kernel_base_cat, 3},
    {"_lgpr_STAN_kernel_base_bin_mask", (DL_FUNC) &_lgpr_STAN_kernel_base_bin_mask, 3},
    {"_lgpr_STAN_kernel_const", (DL_FUNC) &_lgpr_STAN_kernel_const, 5},
    {"_lgpr_STAN_kernel_const_all", (DL_FUNC) &_lgpr_STAN_kernel_const_all, 9},
    {"_lgpr_STAN_kernel_base_var_mask", (DL_FUNC) &_lgpr_STAN_kernel_base_var_mask, 5},
    {"_lgpr_STAN_kernel_all", (DL_FUNC) &_lgpr_STAN_kernel_all, 16},
    {"_lgpr_STAN_gp_posterior_helper", (DL_FUNC) &_lgpr_STAN_gp_posterior_helper, 4},
    {"_lgpr_STAN_gp_posterior", (DL_FUNC) &_lgpr_STAN_gp_posterior, 5},
    {"_lgpr_STAN_log_prior", (DL_FUNC) &_lgpr_STAN_log_prior, 4},
    {"_rcpp_module_boot_stan_fit4lgp_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4lgp_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_lgpr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
