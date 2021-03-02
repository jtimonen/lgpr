  // Binary option switches
  int<lower=0, upper=1> is_verbose;

  // Dimensions
  int<lower=0> num_obs;         // number of observations
  int<lower=0> num_cov_cont;    // number of continuous covariates
  int<lower=0> num_cov_cat;     // number of categorical covariates
  int<lower=1> num_comps;       // number of additive components
  int<lower=0> num_ell;         // number of lengthscale parameters
  int<lower=0> num_ns;          // number of nonstationary components
  int<lower=0> num_heter;       // number of heterogeneous components
  int<lower=0> num_uncrt;       // number of uncertain continuous covariates
  int<lower=0> num_bt;          // number of beta and/or teff params

  /* 
  Each additive function component can be related to at most one continuous and
  one categorical covariate. Properties of the components are specified by the 
  "columns" of the integer array <components>. The "rows" are
    - [,1]: component type
    - [,2]: kernel type
    - [,3]: (currently unused row)
    - [,4]: is the effect magnitude heterogeneous?
    - [,5]: should input warping be applied to the continuous covariate first?
    - [,6]: should a variance mask be applied?
    - [,7]: is there uncertainty in the continuous covariate?
    - [,8]: index of the categorical covariate in <x_cat>, <x_cat_num_levels>
    - [,9]: index of the continuous covariate in <x_cont>, <x_cont_mask>
    
  NOTES: Options [,6] and [,7] only have an effect if option [,5] is 1.
  Possible types for option [,1] are
    - type 0 = component with a single categorical covariate
      * kernel 0 = zero-sum kernel
      * kernel 1 = categorical kernel
    - type 1 = component with a single continuous covariate
      * kernel 0 = [exp. quadratic] kernel
    - type 2 = interaction of a categorical and a continuous covariate
      * kernel 0 = zero-sum kernel * [exp. quadratic]
      * kernel 1 = categorical kernel * [exp. quadratic]
  */
  int<lower=0> components[num_comps, 9];

  // Covariates
  vector[num_obs] x_cont[num_cov_cont];
  vector[num_obs] x_cont_unnorm[num_cov_cont];
  int x_cont_mask[num_cov_cont, num_obs];
  int x_cat[num_cov_cat, num_obs];
  int<lower=0> x_cat_num_levels[num_cov_cat];
  int<lower=1, upper=num_bt+1> idx_expand[num_obs]; // expands beta and t_eff.

  // Priors
  int<lower=0> prior_alpha[num_comps, 2];  // {prior_type, transform}
  int<lower=0> prior_ell[num_ell, 2];      // {prior_type, transform}
  int<lower=0> prior_wrp[num_ns, 2];       // {prior_type, transform}
  int<lower=0> prior_teff[num_uncrt>0, 2]; // {prior type, is_backwards}
  real hyper_alpha[num_comps, 3];
  real hyper_ell[num_ell, 3];
  real hyper_wrp[num_ns, 3];
  real hyper_teff[num_uncrt>0, 3];
  real hyper_beta[num_heter>0, 2];

  // Observed effect times and uncertainty bounds for each case subject
  vector[num_bt] teff_zero[num_uncrt>0];
  vector[num_bt] teff_lb[num_uncrt>0];
  vector[num_bt] teff_ub[num_uncrt>0];

  // Misc
  real delta; // jitter to ensure pos. def. kernel matrices
  real vm_params[num_ns, 2]; // variance mask parameters (nonstat comps)
