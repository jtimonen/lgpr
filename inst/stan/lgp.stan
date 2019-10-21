// *lgp.stan*
// This is the main Stan model file of the 'lgpr' package
// Author: Juho Timonen

/*
    lgpr is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    lgpr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with lgpr.  If not, see <http://www.gnu.org/licenses/>.
*/


// -------------------- 1. FUNCTIONS -------------------- //

functions{

  // LOG PRIOR TO BE ADDED TO TARGET
  real STANFUNC_log_prior(real x, int[ ] types, real[ ] hp){
    real lp;
    real a = hp[1]; // prior hyperparameter 1
    real b = hp[2]; // prior hyperparameter 2
    real c = hp[3]; // prior hyperparameter 3 (currently not used)
    real theta;

    // Possible transform and log of its absolute derivative
    if (types[2]==0){
      lp = 0;
      theta = x;
    }else if (types[2]==1){
      lp = log(fabs(2*x));
      theta = square(x);
    }else{
      reject("invalid value of types[2]!")
    }

    // Value of pdf
    if (types[1]==1){
      // do nothing
    }else if (types[1]==2){
      lp += normal_lpdf(theta|a,b);
    }else if (types[1]==3){
      lp += student_t_lpdf(theta|a,0,b);
    }else if (types[1]==4){
      lp += gamma_lpdf(theta|a,b);
    }else if (types[1]==5){
      lp += inv_gamma_lpdf(theta|a,b);
    }else if (types[1]==6){
      lp += lognormal_lpdf(theta|a,b);
    }else{
      reject("types[1] must be an integer between 1 and 6; found =", types[1]);
    }

    return(lp);
  }

  // INPUT WARP
  vector STANFUNC_warp_input(vector x, real a, real b, real c){
    int n_tot = num_elements(x);
    vector[n_tot] w = 2*c*(-0.5 + rep_vector(1, n_tot)./(1+exp(-a*(x-b))));
    return(w);
  }

  // VARIANCE MASK
  vector STANFUNC_var_mask(vector x, real a){
    int n_tot = num_elements(x);
    vector[n_tot] s = rep_vector(1, n_tot)./(1+exp(-a*x));
    return(s);
  }

  // COMPUTE X_TILDE
  vector STANFUNC_get_x_tilde(vector x_disAge, vector T_onset, vector T_observed, int[,] mapping, int[] lengths){
    int n_tot = num_elements(x_disAge);
    int N_cases = num_elements(lengths);
    vector[n_tot] x_tilde = rep_vector(0.0, n_tot);
    for(k in 1:N_cases){
      int inds[lengths[k]] = mapping[k, 1:lengths[k]];
      x_tilde[inds] = x_disAge[inds] + T_observed[k] - T_onset[k];
    }
    return(x_tilde);
  }

  // CATEGORICAL KERNEL
  matrix STANFUNC_K_cat(vector x1, vector x2){
    int n1 = num_elements(x1);
    int n2 = num_elements(x2);
    matrix[n1,n2] K;
    for(i in 1:n1){
      for(j in 1:n2){
        K[i,j] = x1[i]==x2[j];
      }
    }
    return(K);
  }

  // BINARY (MASK) KERNEL
  matrix STANFUNC_K_bin_int(int[] x1, int[] x2, int c){
    int n1 = num_elements(x1);
    int n2 = num_elements(x2);
    matrix[n1,n2] K;
    for(i in 1:n1){
      for(j in 1:n2){
        K[i,j] = (x1[i]==c)*(x2[j]==c);
      }
    }
    return(K);
  }

  // Binary (mask) kernel but with real inputs
  matrix STANFUNC_K_bin_real(vector x1, vector x2, int c){
    int n1 = num_elements(x1);
    int n2 = num_elements(x2);
    matrix[n1,n2] K;
    for(i in 1:n1){
      for(j in 1:n2){
        K[i,j] = (x1[i]==c)*(x2[j]==c);
      }
    }
    return(K);
  }

  // Multiplier matrix to enable variance masking
  matrix STANFUNC_K_var_mask(vector x1_tilde, vector x2_tilde, real stp, real[] vm_params){
    int n1 = num_elements(x1_tilde);
    int n2 = num_elements(x2_tilde);
    real a = stp * vm_params[2];
    real h = vm_params[1];
    real r = 1/a*log(h/(1-h));
    vector[n1] s1 = STANFUNC_var_mask(x1_tilde - r, a);
    vector[n2] s2 = STANFUNC_var_mask(x2_tilde - r, a);
    matrix[n1,n2] S1 = rep_matrix(s1, n2);
    matrix[n1,n2] S2 = transpose(rep_matrix(s2, n1));
    matrix[n1,n2] MASK = S1 .* S2;
    return(MASK);
  }

  // Multiplier matrix to enable heterogeneous diseaseAge effect
  matrix STANFUNC_K_beta_symmetric(vector beta, int[] row_to_caseID){
    int n = num_elements(row_to_caseID);
    int i_caseID = 0;
    int j_caseID = 0;
    real tmp;
    matrix[n,n] BETA;
    for(i in 1:(n-1)){
      i_caseID = row_to_caseID[i];
      for(j in (i+1):n){
        j_caseID = row_to_caseID[j];
        if(i_caseID*j_caseID > 0){
          tmp = sqrt(beta[i_caseID]*beta[j_caseID]);
        }else{
          tmp = 0;
        }
        BETA[i,j] = tmp;
        BETA[j,i] = tmp;

      }
      if(i_caseID > 0){
        BETA[i,i] = beta[i_caseID];
      }else{
        BETA[i,i] = 0;
      }
    }
    i_caseID = row_to_caseID[n];
    if(i_caseID > 0){
      BETA[n,n] = beta[i_caseID];
    }else{
      BETA[n,n] = 0;
    }
    return(BETA);
  }

    // Multiplier matrix to enable heterogeneous diseaseAge effect
  matrix STANFUNC_K_beta(vector beta, int[] row_to_caseID_1, int[] row_to_caseID_2){
    int n1 = num_elements(row_to_caseID_1);
    int n2 = num_elements(row_to_caseID_2);
    int i_caseID = 0;
    int j_caseID = 0;
    real tmp;
    matrix[n1,n2] BETA;
    for(i in 1:n1){
      i_caseID = row_to_caseID_1[i];
      for(j in 1:n2){
        j_caseID = row_to_caseID_2[j];
        if(i_caseID*j_caseID > 0){
          BETA[i,j] = sqrt(beta[i_caseID]*beta[j_caseID]);
        }else{
          BETA[i,j] = 0;
        }
      }
    }
    return(BETA);
  }


  // COMPUTE FIXED KERNEL MATRICES (do not depend on parameters)
  matrix[] STANFUNC_compute_fixed_kernel_matrices(vector[] X1, vector[] X2, int[] X1_nn, int[] X2_nn, int[] D, int cat_interact_kernel){
    int n1 = num_elements(X1[1]);
    int n2 = num_elements(X2[1]);
    int nf = 1 + D[3] + D[5] + D[6];
    matrix[n1,n2] KF[nf];
    KF[1] = STANFUNC_K_cat(X1[1], X2[1]);
    for(j in 1:D[3]){
      KF[1+j] = STANFUNC_K_bin_int(X1_nn, X2_nn, 1);
    }
    for(j in 1:D[5]){
      int ix = 2 + D[3] + D[4] + j;
      if(cat_interact_kernel == 1){
        KF[1+D[3]+j] = STANFUNC_K_cat(X1[ix], X2[ix]);
      }else{
        KF[1+D[3]+j] = STANFUNC_K_bin_real(X1[ix], X2[ix], 1);
      }
    }
    for(j in 1:D[6]){
      int ix = 2+D[3]+D[4]+D[5]+j;
      KF[1+D[3]+D[5]+j] = STANFUNC_K_cat(X1[ix], X2[ix]);
    }
    return(KF);
  }


  // COMPUTE ALL KERNEL MATRICES
  matrix[] STANFUNC_compute_kernel_matrices(vector[] X1, vector[] X2, int[,] caseID_to_rows_1, int[,] caseID_to_rows_2, int[] row_to_caseID_1, int[] row_to_caseID_2, int[] caseID_nrows_1, int[] caseID_nrows_2, matrix[] KF, vector[] T_onset, vector T_observed, int[] D, int UNCRT, int HMGNS, int USE_VAR_MASK, real[] vm_params, real[] alpha_idAge, real[] alpha_sharedAge, real[] alpha_diseaseAge, real[] alpha_continuous, real[] alpha_categAge, real[] alpha_categOffset, real[] lengthscale_idAge, real[] lengthscale_sharedAge, real[] lengthscale_diseaseAge, real[] lengthscale_continuous, real[] lengthscale_categAge, real[] warp_steepness, vector[] beta){
    int n1 = num_elements(X1[1]);
    int n2 = num_elements(X2[1]);
    real x1_age[n1] = to_array_1d(X1[2]);   // age covariate as an array
    real x2_age[n2] = to_array_1d(X2[2]);   // age covariate as an array
    int sum_D = sum(D);
    matrix[n1,n2] KX[sum_D];
    int r = 0;
    if(D[1]==1){
      real alp = alpha_idAge[1];
      real ell = lengthscale_idAge[1];
      r += 1;
      KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell) .* KF[1];
    }
    if(D[2]==1){
      real alp = alpha_sharedAge[1];
      real ell = lengthscale_sharedAge[1];
      r += 1;
      KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell);
    }
    if(D[3]==1){
      real alp = alpha_diseaseAge[1];
      real ell = lengthscale_diseaseAge[1];
      real stp = warp_steepness[1];
      vector[n1] x1_tilde;
      vector[n2] x2_tilde;
      real w1[n1];
      real w2[n2];
      r += 1;

      // Handle diseaseAge uncertainty
      if(UNCRT==0){
        x1_tilde = X1[3];
        x2_tilde = X2[3];
      }else{
        x1_tilde = STANFUNC_get_x_tilde(X1[3], T_onset[1], T_observed, caseID_to_rows_1, caseID_nrows_1);
        x2_tilde = STANFUNC_get_x_tilde(X2[3], T_onset[1], T_observed, caseID_to_rows_2, caseID_nrows_2);
      }
      w1 = to_array_1d(STANFUNC_warp_input(x1_tilde, stp, 0.0, 1.0));
      w2 = to_array_1d(STANFUNC_warp_input(x2_tilde, stp, 0.0, 1.0));

      // Create disease effect kernel
      KX[r] = KF[2] .* cov_exp_quad(w1, w2, alp, ell);
      if(HMGNS==0){
        KX[r] = STANFUNC_K_beta(beta[1], row_to_caseID_1, row_to_caseID_2) .* KX[r];
      }
      if(USE_VAR_MASK==1){
        KX[r] = STANFUNC_K_var_mask(x1_tilde, x2_tilde, stp, vm_params) .* KX[r];
      }
    }
    for(j in 1:D[4]){
      real alp = alpha_continuous[j];
      real ell = lengthscale_continuous[j];
      int ix  = 2 + D[3] + j;
      r += 1;
      KX[r] = cov_exp_quad(to_array_1d(X1[ix]), to_array_1d(X2[ix]), alp, ell);
    }
    for(j in 1:D[5]){
      real alp = alpha_categAge[j];
      real ell = lengthscale_categAge[j];
      int ikf = 1 + D[3] + j;
      r += 1;
      KX[r] = cov_exp_quad(x1_age, x2_age, alp, ell) .* KF[ikf];
    }
    for(j in 1:D[6]){
      int ikf = 1 + D[3] + D[5] + j;
      real alp = alpha_categOffset[j];
      r += 1;
      KX[r] = square(alp) * KF[ikf];
    }
    return(KX);
  }

}



// -------------------- 2. DATA  -------------------- //

data {

  int<lower=1> N_tot;        // total number of individuals
  int<lower=0> N_cases;      // number of "diseased" individuals
  int<lower=2> d;            // number of covariates in the data (id and age required)
  int<lower=1> n;            // number of observations
  int<lower=1> N_trials[n];  // numbers of trials (set to all ones for bernoulli model)

  // Modeled data
  vector[n] X[d];           // covariates, X[j] is the jth covariate
  int       X_id[n];        // the id covariate as an array of integers
  int       X_notnan[n];    // X_notnan[i] tells if X_diseaseAge[i] is originally NaN
  vector[n] y;              // the response variable (as a vector of reals)
  int       y_int[n];       // the response variable (as an array of integers)

  // Observation model
  int<lower=0,upper=4> LH;

  // D is an array of six integers, so that
  //   D[1] = binary value indicating if id*age is a predictor
  //   D[2] = binary value indicating if shared age is a predictor
  //   D[3] = binary value indicating if diseaseAge is a predictor
  //   D[4] = number of other continuous covariates
  //   D[5] = number of discrete covariates that interact with age
  //   D[6] = number of discrete covariates that only have an offset effect
  int<lower=0> D[6];

  // Modeling option switches
  int<lower=0,upper=1> UNCRT;  // tells if the diseaseAge measurements are uncertain
  int<lower=0,upper=1> HMGNS;  // tells if the diseaseAge effect is homogenous

  // Prior types and transforms
  int t_ID[D[1],4];         // for id*age component
  int t_A[D[2],4];          // for shared age component
  int t_D[D[3],6];          // for disease age component
  int t_CNT[D[4],4];        // for components with continuous covariates
  int t_CAT[D[5],4];        // for components with discrete covariates
  int t_OFS[D[6],2];        // for offset components
  int t_SIG[2];             // for Gaussian noise std
  int t_PHI[2];             // for precision parameter phi
  int t_ONS[N_cases,2];     // for onset, if uncertain

  // Hyperparameters of the priors and scaling factors,
  real p_ID[D[1],6];         // for id*age component
  real p_A[D[2],6];          // for shared age component
  real p_D[D[3],9];          // for disease age component
  real p_CNT[D[4],6];        // for components with continuous covariates
  real p_CAT[D[5],6];        // for components with discrete covariates
  real p_OFS[D[6],3];        // for offset components
  real p_SIG[3];             // for Gaussian noise std
  real p_PHI[3];             // for precision parameter phi
  real p_BET[2];             // hyperparameters of the beta prior of BETA
  real p_ONS[N_cases, 3];    // for onset, if uncertain

  // Inputs related to mapping from row index to case index and back
  int<lower=0,upper=n>        M_max;
  int<lower=0>                caseID_to_rows[N_cases, M_max];
  int<lower=0,upper=N_cases>  row_to_caseID[n];
  int<lower=0,upper=M_max>    caseID_nrows[N_cases];

  // Inputs related to uncertain disease onset
  vector[N_cases] T_observed;      // observed disease effect times
  vector[N_cases] L_ons[UNCRT];    // lower bounds for sampled disease effect times
  vector[N_cases] U_ons[UNCRT];    // upper bounds for sampled disease effect times
  int<lower=0,upper=1> backwards;  // is the prior for onset "backwards"
  int<lower=0,upper=1> relative;   // is the prior for effect time relative to observed one

  // Other
  real DELTA;       // jitter to ensure pos. def. kernel matrices
  real C_hat;       // C_hat parameter for poisson and NB models
  int HS[6];        // (currently not used)
  int F_is_sampled; // should the function values be sampled?
                    // (must be 1 for non-Gaussian lh)

  // Kernel types that have options
  int<lower=0,upper=1> USE_VAR_MASK;
  real vm_params[2];
  int<lower=0,upper=1> cat_interact_kernel; // {1 = categorical, 0 = binary} kernel

}



// -------------------- 3. TRANSFORMED DATA  -------------------- //

transformed data{
  int nf = 1 + D[3] + D[5] + D[6];     // number of fixed kernel matrices
  int sum_D = sum(D);                  // total number of covariates
  matrix[n,n] KF[nf] = STANFUNC_compute_fixed_kernel_matrices(X, X, X_notnan, X_notnan, D, cat_interact_kernel);
  vector[n] mu = rep_vector(C_hat, n); // GP mean


  // Print some info (mostly for debugging)
  print(" ")
  print("* Observation model = ", LH);
  print("* Number of data points = ", n);
  print("* Number of model components = ", sum_D);
  print("* Number of individuals = ", N_tot);
  print("* Additional model info:")
  if(LH==2 || LH==3){
    print("  - C_hat = ", C_hat);
  }
  print("  - D = ", D);
  print("  - F_is_sampled = ", F_is_sampled)
  print("  - cat_interact_kernel = ", cat_interact_kernel)
  if(D[3]==1){
    print("* Disease modeling info: ");
    print("  - Number of cases = ", N_cases);
    print("  - UNCRT = ", UNCRT);
    print("  - HMGNS = ", HMGNS);
    print("  - USE_VAR_MASK = ", USE_VAR_MASK);
    if(USE_VAR_MASK==1){
      print("      o vm_params = ", vm_params);
    }
  }
  print(" ")

}



// -------------------- 4. PARAMETERS  -------------------- //

parameters {

  // Magnitude params
  real<lower=0> alpha_idAge[D[1]];
  real<lower=0> alpha_sharedAge[D[2]];
  real<lower=0> alpha_diseaseAge[D[3]];
  real<lower=0> alpha_continuous[D[4]];
  real<lower=0> alpha_categAge[D[5]];
  real<lower=0> alpha_categOffset[D[6]];

  // Lengthscale parameters
  real<lower=0> lengthscale_idAge[D[1]];
  real<lower=0> lengthscale_sharedAge[D[2]];
  real<lower=0> lengthscale_diseaseAge[D[3]];
  real<lower=0> lengthscale_continuous[D[4]];
  real<lower=0> lengthscale_categAge[D[5]];

  // Miscellaneous
  real<lower=0> warp_steepness[D[3]];       // steepness of input warping
  real<lower=0> sigma_n[LH==1 || LH==0];    // noise std for Gaussian likelihood
  vector[n] ETA[F_is_sampled, sum_D];       // isotropic versions of F
  real<lower=0> phi[LH==3 || LH==0];        // overdispersion parameter for NB likelihood

  // Parameters related to diseased individuals
  vector<lower=0,upper=1>[N_cases] beta[HMGNS==0];  // individual-specific magnitude
  vector<lower=0,upper=1>[N_cases] T_raw[UNCRT==1];

}



// -------------------- 5. TRANSFORMED PARAMETERS  -------------------- //

transformed parameters {
  vector[N_cases] T_onset[UNCRT];
  vector[n] F[F_is_sampled, sum_D];
  if(UNCRT){
    T_onset[1] = L_ons[1] + (U_ons[1] - L_ons[1]) .* T_raw[1];
  }
  if(F_is_sampled){
    matrix[n,n] KX[sum_D] = STANFUNC_compute_kernel_matrices(X, X, caseID_to_rows, caseID_to_rows, row_to_caseID, row_to_caseID, caseID_nrows, caseID_nrows, KF, T_onset, T_observed, D, UNCRT, HMGNS, USE_VAR_MASK, vm_params, alpha_idAge, alpha_sharedAge,  alpha_diseaseAge, alpha_continuous, alpha_categAge, alpha_categOffset, lengthscale_idAge, lengthscale_sharedAge, lengthscale_diseaseAge, lengthscale_continuous, lengthscale_categAge, warp_steepness, beta);
    for(r in 1:sum_D){
      matrix[n,n] EYE = diag_matrix(rep_vector(DELTA, n));
      matrix[n,n] Lxr = cholesky_decompose(KX[r] + EYE);
      F[1,r,] = Lxr*ETA[1,r,];
    }
  }
}



// -------------------- 6. MODEL  -------------------- //

model {

  // 6.1 PRIORS

  // Kernel hyperparameters for idAge component
  if(D[1]==1){
    target += STANFUNC_log_prior(alpha_idAge[1],       t_ID[1,1:2], p_ID[1,1:3]);
    target += STANFUNC_log_prior(lengthscale_idAge[1], t_ID[1,3:4], p_ID[1,4:6]);
  }

  // Kernel hyperparameters for sharedAge component
  if(D[2]==1){
    target += STANFUNC_log_prior(alpha_sharedAge[1],       t_A[1,1:2], p_A[1,1:3]);
    target += STANFUNC_log_prior(lengthscale_sharedAge[1], t_A[1,3:4], p_A[1,4:6]);
  }

  // Kernel hyperparameters for diseaseAge component
  if(D[3]==1){
    target += STANFUNC_log_prior(alpha_diseaseAge[1],       t_D[1,1:2], p_D[1,1:3]);
    target += STANFUNC_log_prior(lengthscale_diseaseAge[1], t_D[1,3:4], p_D[1,4:6]);
    target += STANFUNC_log_prior(warp_steepness[1],         t_D[1,5:6], p_D[1,7:9]);
  }

  // Kernel hyperparameters for other continuous components
  for(j in 1:D[4]){
    target += STANFUNC_log_prior(alpha_continuous[j],       t_CNT[j,1:2], p_CNT[j,1:3]);
    target += STANFUNC_log_prior(lengthscale_continuous[j], t_CNT[j,3:4], p_CNT[j,4:6]);
  }

  // Kernel hyperparameters for categorical * age components
  for(j in 1:D[5]){
    target += STANFUNC_log_prior(alpha_categAge[j],       t_CAT[j,1:2], p_CAT[j,1:3]);
    target += STANFUNC_log_prior(lengthscale_categAge[j], t_CAT[j,3:4], p_CAT[j,4:6]);
  }

  // Kernel hyperparameters for categorical offset components
  for(j in 1:D[6]){
    target += STANFUNC_log_prior(alpha_categOffset[j], t_OFS[j,1:2], p_OFS[j,1:3]);
  }

  // Noise level parameters
  if(LH==1){
    target += STANFUNC_log_prior(sigma_n[1], t_SIG[1:2], p_SIG[1:3]);
  }else if(LH==3){
    target += STANFUNC_log_prior(phi[1], t_PHI[1:2], p_PHI[1:3]);
  }else if(LH==0){
    target += STANFUNC_log_prior(sigma_n[1], t_SIG[1:2], p_SIG[1:3]);
    target += STANFUNC_log_prior(phi[1],     t_PHI[1:2], p_PHI[1:3]);
  }

  // DiseaseAge uncertainty prior
  if(UNCRT){
    real tx;
    for(k in 1:N_cases){
      if(relative==1){
        tx = - T_observed[k] + T_onset[1,k];
      }else if(backwards==1){
        tx = T_observed[k] - T_onset[1,k];
      }else{
        tx = T_onset[1,k];
      }
      target += STANFUNC_log_prior(tx, t_ONS[k,1:2], p_ONS[k,1:3]);
    }
  }

  // Heterogeneous disease effect parameters
  if(HMGNS==0){
    beta[1] ~ beta(p_BET[1], p_BET[2]);
  }

  // Isotropic normals
  if(F_is_sampled){
    for(j in 1:sum_D){
      ETA[1,j] ~ normal(0,1);
    }
  }

  // 6.2 OBSERVATION MODEL
  if(LH==0){

  }else{

    if(F_is_sampled){

      //  F IS SAMPLED
      // Compute f
      vector[n] F_sum = rep_vector(0, n);
      for (i in 1:n){
        F_sum[i] += sum(F[1,,i]);
      }

      // Compute likelihood
      if(LH==1){
        // 1. Gaussian observation model
        real SIGMA[n] = to_array_1d(rep_vector(sigma_n[1], n)); // means
        real MU[n] = to_array_1d(F_sum[1:n]);                   // stds
        target += normal_lpdf(y | MU, SIGMA);
      }else if(LH==2){
        // 2. Poisson observation model
        real LOG_MU[n] = to_array_1d(F_sum[1:n] + C_hat); // means (rate parameters) on log-scale
        target += poisson_log_lpmf(y_int | LOG_MU);
      }else if(LH==3){
        // 3. Negative binomial observation model
        real LOG_MU[n] = to_array_1d(F_sum[1:n] + C_hat); // means on log-scale
        real PHI[n] = to_array_1d(rep_vector(phi[1], n)); // inverse dispersion parameters
        target += neg_binomial_2_log_lpmf(y_int | LOG_MU, PHI);
      }else if(LH==4){
        // 4. Bernoulli or binomial observation model
        real LOGIT_P[n] = to_array_1d(F_sum[1:n]); // probabilities of success on log-scale
        target += binomial_logit_lpmf(y_int | N_trials, LOGIT_P);
      }
      else{
        reject("Unknown observation model!")
      }

    }else{
      // F NOT SAMPLED
      matrix[n,n] Ky;
      matrix[n,n] Ly;
      matrix[n,n] Kx = diag_matrix(rep_vector(DELTA, n));
      matrix[n,n] KX[sum_D] = STANFUNC_compute_kernel_matrices(X, X, caseID_to_rows, caseID_to_rows, row_to_caseID, row_to_caseID, caseID_nrows, caseID_nrows, KF, T_onset, T_observed, D, UNCRT, HMGNS, USE_VAR_MASK, vm_params, alpha_idAge, alpha_sharedAge,  alpha_diseaseAge,  alpha_continuous, alpha_categAge, alpha_categOffset, lengthscale_idAge, lengthscale_sharedAge, lengthscale_diseaseAge, lengthscale_continuous, lengthscale_categAge, warp_steepness, beta);
      if(LH!=1){
        reject("Likelihood must be Gaussian if F is not sampled!")
      }
      for(j in 1:sum_D){
        Kx += KX[j];
      }
      Ky = Kx + diag_matrix(rep_vector(square(sigma_n[1]), n));
      Ly = cholesky_decompose(Ky);
      y ~ multi_normal_cholesky(mu, Ly);
    }
  }
}



// -------------------- 7. GENERATED QUANTITIES  -------------------- //

generated quantities {
   vector[n] F_mean_cmp[1 - F_is_sampled, sum_D];
   vector[n] F_var_cmp[ 1 - F_is_sampled, sum_D];
   vector[n] F_mean_tot[1 - F_is_sampled];
   vector[n] F_var_tot[ 1 - F_is_sampled];

   if(F_is_sampled==0){
     matrix[n,n] A;
     vector[n] v;
     matrix[n,n] Ky;
     matrix[n,n] Ly;
     matrix[n,n] Kx = diag_matrix(rep_vector(DELTA, n));
     matrix[n,n] KX[sum_D] = STANFUNC_compute_kernel_matrices(X, X, caseID_to_rows, caseID_to_rows, row_to_caseID, row_to_caseID, caseID_nrows, caseID_nrows, KF, T_onset, T_observed, D, UNCRT, HMGNS, USE_VAR_MASK, vm_params, alpha_idAge, alpha_sharedAge,  alpha_diseaseAge,  alpha_continuous, alpha_categAge, alpha_categOffset, lengthscale_idAge, lengthscale_sharedAge, lengthscale_diseaseAge, lengthscale_continuous, lengthscale_categAge, warp_steepness, beta);
     for(j in 1:sum_D){
       Kx += KX[j];
     }
     Ky = Kx + diag_matrix(rep_vector(square(sigma_n[1]), n));
     Ly = cholesky_decompose(Ky);
     v = mdivide_left_tri_low(Ly, y);
     for(j in 1:sum_D){
       A  = mdivide_left_tri_low(Ly, transpose(KX[j]));
       F_mean_cmp[1,j] = transpose(A)*v;
       F_var_cmp[1,j] = diagonal(KX[j] - crossprod(A));
     }
     A = mdivide_left_tri_low(Ly, transpose(Kx));
     F_mean_tot[1] = transpose(A)*v;
     F_var_tot[1] = diagonal(Kx - crossprod(A));
   }
}

