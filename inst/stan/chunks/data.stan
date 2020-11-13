/* 
  Binary option switches
  - is verbose mode used?
  - are function values be sampled?
  - should likelihood evaluation be skipped?
*/
int<lower=0, upper=1> is_verbose;
int<lower=0, upper=1> is_f_sampled;
int<lower=0, upper=1> is_likelihood_skipped;

// Dimensions
int<lower=0> num_obs;            // number of observations
int<lower=0> num_cov_cont;       // number of continuous covariates
int<lower=0> num_cov_cat;        // number of categorical covariates
int<lower=1> num_comps;          // number of additive components
int<lower=0> num_ell;            // number of ell parameters
int<lower=0> num_ns;             // number of nonstationary components
int<lower=0> num_heter;          // number of heterogeneous components
int<lower=0> num_uncrt;          // number of uncertain continuous covariates
int<lower=0> num_bt;             // number of beta and/or teff params

/*
  Observation model
  - 1 = Gaussian
  - 2 = Poisson
  - 3 = Negative Binomial
  - 4 = Binomial
  - 5 = Beta-Binomial
*/
int<lower=1,upper=5> obs_model;

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

// Response variable (vector of reals)
vector[num_obs] y_cont[obs_model==1];

// Response variable (array of integers)
int<lower=0> y_disc[obs_model>1, num_obs];

// Covariates
vector[num_obs] x_cont[num_cov_cont];
vector[num_obs] x_cont_unnorm[num_cov_cont];
int x_cont_mask[num_cov_cont, num_obs];
int x_cat[num_cov_cat, num_obs];

// Number of trials (binomial or bernoulli model)
int<lower=1> y_num_trials[obs_model>3, num_obs];

// Number of levels for each categorical covariate
int<lower=0> x_cat_num_levels[num_cov_cat];

// Inputs related to expanding beta and t_effect
int<lower=1, upper=num_bt+1> idx_expand[num_obs];

/* 
  Prior types and transforms for positive parameters
  - [,1]: prior types
  - [,2]: prior transforms
*/
int<lower=0> prior_alpha[num_comps, 2];
int<lower=0> prior_ell[num_ell, 2];
int<lower=0> prior_wrp[num_ns, 2];
int<lower=0> prior_sigma[obs_model==1, 2];
int<lower=0> prior_phi[obs_model==3, 2];

/*
  Prior types for effect time uncertainty
  - [,1]: prior type
  - [,2]: is prior "backwards"?
*/
int<lower=0> prior_teff[num_uncrt>0, 2];

// Hyperparameters of the priors for positive parameters
real hyper_alpha[num_comps, 3];
real hyper_ell[num_ell, 3];
real hyper_wrp[num_ns, 3];
real hyper_sigma[obs_model==1, 3];
real hyper_phi[obs_model==3, 3];
real hyper_teff[num_uncrt>0, 3];

// Hyperparameters of the beta priors for parameters on [0, 1]
real hyper_gamma[obs_model==5, 2];
real hyper_beta[num_heter>0, 2];

// Observed effect times and uncertainty bounds for each case subject
vector[num_bt] teff_zero[num_uncrt>0];
vector[num_bt] teff_lb[num_uncrt>0];
vector[num_bt] teff_ub[num_uncrt>0];

// Misc
vector[num_obs] c_hat; // "GP mean vector"
real delta; // jitter to ensure pos. def. kernel matrices
real vm_params[num_ns, 2]; // variance mask parameters

// Basis function approximation
int<lower=0> num_basisfun; // 0 = no approximation used
real<lower=0> width_basisfun;
