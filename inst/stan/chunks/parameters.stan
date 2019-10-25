
// Magnitude params
real<lower=0> alpha_idAge[D[1]];
real<lower=0> alpha_sharedAge[D[2]];
real<lower=0> alpha_diseaseAge[D[3]];
real<lower=0> alpha_continuous[D[4]];
real<lower=0> alpha_categAge[D[5]];
real<lower=0> alpha_categOffset[D[6]];

// Lengthscale parameters
real<lower=0> ell_idAge[D[1]];
real<lower=0> ell_sharedAge[D[2]];
real<lower=0> ell_diseaseAge[D[3]];
real<lower=0> ell_continuous[D[4]];
real<lower=0> ell_categAge[D[5]];

// Miscellaneous
real<lower=0> warp_steepness[D[3]];     // steepness of input warping
real<lower=0> sigma_n[LH==1 || LH==0];  // noise std for Gaussian likelihood
vector[n] ETA[F_IS_SAMPLED, sum_D];     // isotropic versions of F
real<lower=0> phi[LH==3 || LH==0];      // parameter for NB likelihood

// Parameters related to diseased individuals
vector<lower=0,upper=1>[N_cases] beta[HMGNS==0]; 
vector<lower=0,upper=1>[N_cases] T_raw[UNCRT==1];
