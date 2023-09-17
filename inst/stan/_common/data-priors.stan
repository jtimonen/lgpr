  // Priors types and hyperparameters
  array[num_comps, 2] int<lower=0> prior_alpha;  // {prior_type, transform}
  array[num_ell, 2] int<lower=0> prior_ell;      // {prior_type, transform}
  array[num_ns, 2] int<lower=0> prior_wrp;       // {prior_type, transform}
  array[num_uncrt>0, 2] int<lower=0> prior_teff; // {prior type, is_backwards}
  array[num_comps, 3] real hyper_alpha;
  array[num_ell, 3] real hyper_ell;
  array[num_ns, 3] real hyper_wrp;
  array[num_uncrt>0, 3] real hyper_teff;
  array[num_heter>0, 2] real hyper_beta;
