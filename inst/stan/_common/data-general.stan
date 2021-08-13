  // Binary option switches
  int<lower=0, upper=1> is_verbose;
  int<lower=0, upper=1> is_likelihood_skipped;
  
  // Dimensions
  int<lower=0> N;           // number of observations
  int<lower=1> J;           // number of components
  int<lower=0> num_X;       // number of continuous covariates
  int<lower=0> num_Z;       // number of categorical covariates
  int<lower=0> num_ell;     // number of lengthscale parameters
  int<lower=0> num_wrp;     // number of input warping parameters
  int<lower=0> num_beta;    // number of beta params
  int<lower=0> num_teff;    // number of effect time params
  int<lower=0, upper=num_X> idx_unc; // index of uncertain covariate (0 = none)
  int<lower=0> num_het;     // number of heterogeneous components

  /*
  Each additive function component can be related to at most one continuous and
  one categorical covariate. Properties of the components are specified by the 
  "columns" of the integer array <components>. The "rows" are
    - [,1]: type of categorical kernel (0 = none, 1 = CAT, 2 = ZS)
    - [,2]: type of continuous kernel (0 = none, 1 = EQ)
    - [,3]: index of the categorical covariate in <Z>, <Z_num_levels>
    - [,4]: index of the continuous covariate in <X>, <X_mask>
    - [,5]: input warping (0 = no, 1 = yes, 2 = yes with variance mask)
    - [,6]: heterogeneity (0 = no, 1 = yes)
  */
  int<lower=0> components[J, 6];

  // Observed effect times and uncertainty bounds
  vector[num_teff] teff_zero[idx_unc>0];
  vector[num_teff] teff_lb[idx_unc>0];
  vector[num_teff] teff_ub[idx_unc>0];

  // Misc
  real delta; // jitter to ensure pos. def. kernel matrices
  real vm_params[2]; // variance mask parameters (nonstat comps)
