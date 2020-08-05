// Binary option switches
//  [1]: is verbose mode used?
//  [2]: are function values be sampled?
//  [3]: is the possible disease effect modeled heterogeneously?
//  [4]: is the possible disease effect time uncertain?
//  [5]: is variance masking used for the possible disease component?
//  [6]: should likelihood evaluation be skipped?
//  [7]: should the generated quantities block be skipped?
int<lower=0, upper=1> is_verbose;
int<lower=0, upper=1> is_f_sampled;
int<lower=0, upper=1> is_heter;
int<lower=0, upper=1> is_uncrt;
int<lower=0, upper=1> is_vm_used;
int<lower=0, upper=1> is_likelihood_skipped;
int<lower=0, upper=1> is_generated_skipped;

// Dimensions
int<lower=0> num_subjects;      // total number of subjects
int<lower=0> num_cases;         // number of case subjects
int<lower=0> num_obs;           // number of observations
int<lower=0> num_cov_cont;      // number of continuous covariates
int<lower=0> num_cov_disc;      // number of discrete covariates
int<lower=0> num_comps;         // number of additive components
int<lower=0> num_ell;           // number of ell parameters
int<lower=0> num_dis;           // number of disease components

// Observation model
// - 1 = Gaussian
// - 2 = Poisson
// - 3 = Negative Binomial
// - 4 = Binomial
int<lower=1,upper=4> obs_model;

// Component types are specified by the first two "rows" of the integer array
// <components>, so that on each "column"
//
// - the first number specifies component type
// - the second number specifies kernel type
//
// The possible types are
//
//    type 0 = component with a single categorical covariate
//        kernel 0 = zero-sum kernel
//        kernel 1 = categorical kernel
//        kernel 2 = binary mask kernel (for category 1)
//
//    type 1 = continuous covariate modeled with a stationary kernel
//        kernel 0 = [exp. quadratic] kernel
//
//    type 2 = interaction of categorical and continuous covariate
//        kernel 0 = zero-sum kernel * [exp. quadratic]
//        kernel 1 = categorical kernel * [exp. quadratic]
//        kernel 2 = binary mask kernel (for category 1) * [exp. quadratic] 
//
//    type 3 = disease component (nonstationary)
//        kernel 0 = input warping
//
// Covariates of each component are specified by the last two "rows" of the
// integer array <components>. The third row specifies an index of a discrete
// covariate in <x_disc> the fourth row specifies an index of a continuous
// covariate in <x_cont>. For a disease component (type 3), the discrete
// covariate is "case id" and the continuous covariate is "disease-related age".
int<lower=0> components[4, num_comps];

// Response variable (vector of reals)
vector[num_obs] y_cont[obs_model==1];

// Response variable (array of integers)
int<lower=0> y_disc[obs_model>1, num_obs];

// Covariates
vector[num_obs] x_cont[num_cov_cont]; // continuous covariates
int x_disc[num_cov_disc, num_obs]; // discrete covariates

// Number of trials (binomial or bernoulli model)
int<lower=1> y_num_trials[obs_model==4, num_obs];

// Number of levels for each discrete covariate
int<lower=0> num_levels[num_cov_disc];

// Inputs related to expanding beta and t_effect
int<lower=1, upper=num_cases+1> idx_expand[(is_heter || is_uncrt), num_obs];

// Prior types and transforms for kernel and noise parameters
//   [,1]: prior types
//   [,2]: prior transforms
int<lower=0> prior_alpha[num_comps, 2];
int<lower=0> prior_ell[num_ell, 2];
int<lower=0> prior_wrp[num_dis, 2];
int<lower=0> prior_sigma[obs_model==1, 2];
int<lower=0> prior_phi[obs_model==3, 2];

// Prior types for effect time uncertainty
//   [,1]: prior type
//   [,2]: is prior "backwards"?
//   [,3]: is prior relative to the observed effect time?
int<lower=0> prior_teff[is_uncrt, 3];

// Hyperparameters of the priors
real hyper_alpha[num_comps, 3];
real hyper_ell[num_ell, 3];
real hyper_wrp[num_dis, 3];
real hyper_sigma[obs_model==1, 3];
real hyper_phi[obs_model==3, 3];
real hyper_teff[is_uncrt, 3];
real hyper_beta[is_heter, 2];

// Observed effect times and uncertainty bounds for each case subject
vector[num_cases] teff_obs[is_uncrt];
vector[num_cases] teff_lb[is_uncrt];
vector[num_cases] teff_ub[is_uncrt];

// Misc
vector[num_obs] c_hat; // GP mean vector 
real delta; // jitter to ensure pos. def. kernel matrices
real vm_params[is_vm_used, 2]; // variance mask parameters
