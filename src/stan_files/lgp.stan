// *lgp.stan*
// This is the main Stan model of the 'lgpr' package
// Author: Juho Timonen

#include /chunks/license.stan

functions{
#include /chunks/kernels_base.stan
#include /chunks/define_prior.stan
}

data {

  // Dimensions
  int<lower=1> N_tot;      // total number of individuals
  int<lower=0> N_cases;    // number of "diseased" individuals
  int<lower=2> d;          // number of covariates in the data (id and age required)
  int<lower=1> n;          // number of observations
  int<lower=0> n_test;     // number of test points
  

  // Modeled data
  vector[n + n_test] X[d];          // covariates, X[j] is the jth covariate
  int       X_id[n + n_test];       // the id covariate as an array of integers
  int       X_notnan[n + n_test];   // X_notnan[i] tells if X_diseaseAge[i] is originally NaN
  vector[n] y;                      // the response variable (as a vector of reals)
  int       y_int[n];               // the response variable (as an array of integers)

  // Likelihood type
  int<lower=0,upper=3> LH;          
  
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
  int<lower=0,upper=n+n_test> M_max;
  int<lower=0>                caseID_to_rows[N_cases, M_max];
  int<lower=0,upper=N_cases>  row_to_caseID[n + n_test]; 
  int<lower=0,upper=M_max>    caseID_nrows[N_cases];
  
  // Inputs related to uncertain disease onset
  vector[N_cases] T_observed;      // observed disease onset times
  vector[N_cases] L_ons[UNCRT];    // lower bounds for sampled disease onset
  vector[N_cases] U_ons[UNCRT];    // upper bounds for sampled disease onset
  int<lower=0,upper=1> backwards;  // is the prior for onset "backwards"
  
  // Other
  real DELTA;       // jitter to ensure pos. def. kernel matrices
  real C_hat;       // GP mean (constant value)
  int HS[6];        // (currently not used)
  int F_is_sampled; // should the function values be sampled? 
                    // (must be 1 for non-Gaussian lh)
                    
  // Variance mask
  int<lower=0,upper=1> USE_VAR_MASK;
  real vm_params[2];
  
}

transformed data{
  int n_tot = n + n_test;                  // total number of points
  int nf = 1 + D[3] + D[5] + D[6];         // number of fixed kernel matrices
  int sum_D = sum(D);                      // total number of covariates
  matrix[n_tot,n_tot] KF[nf];              // declare the array of fixed kernel matrices
  real x_age[n_tot] = to_array_1d(X[2]);   // age covariate as an array

  // GP mean
  vector[n_tot] mu = rep_vector(C_hat, n_tot);
    
  // Precompute fixed kernel matrices
  KF[1]  = K_cat(X[1], X[1]);
  for(j in 1:D[3]){
    KF[1+j] = K_bin(X_notnan, X_notnan, 1);
  }
  for(j in 1:D[5]){
    int ix = 2 + D[3] + D[4] + j;
    KF[1+D[3]+j] = K_cat(X[ix], X[ix]);
  }
  for(j in 1:D[6]){
    int ix = 2+D[3]+D[4]+D[5]+j;
    KF[1+D[3]+D[5]+j] = K_cat(X[ix], X[ix]);
  }
  
  // Print some info (mostly for debugging)
  print(" ")
  print("* Likelihood = ", LH);
  print("* Number of data points = ", n);
  print("* Number of test points = ", n_test);
  print("* Number of model components = ", sum_D);
  print("* Number of individuals = ", N_tot);
  print("* Additional model info:")
  if(LH!=1){
    print("  - C_hat = ", C_hat);
  }
  print("  - D = ", D);
  print("  - F_is_sampled = ", F_is_sampled)
  if(D[3]==1){
    print("* Disease modeling info: ");
    print("  - Number of cases = ", N_cases);
    print("  - UNCRT = ", UNCRT);
    print("  - HMGNS = ", HMGNS);
    print("  - USE_VAR_MASK = ", USE_VAR_MASK);
    if(USE_VAR_MASK==1){
      print("      o vm_params = [", vm_params,"]");
    }
  }
  print(" ")
  
  // Input check
  if((n_test) != 0 && (1-F_is_sampled)){
    reject("Number of test points must be zero if F is not sampled!");
  }
}

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
  real<lower=0> warp_steepness[D[3]];               // steepness of input warping
  real<lower=0> sigma_n[LH==1 || LH==0];            // noise std for Gaussian likelihood
  vector[n_tot] ETA[F_is_sampled, sum_D];           // isotropic versions of F
  real<lower=0> phi[LH==3 || LH==0];                // overdispersion parameter for NB likelihood

  // Parameters related to diseased individuals
  vector<lower=0,upper=1>[N_cases] beta[HMGNS==0];  // individual-specific magnitude
  vector<lower=0,upper=1>[N_cases] T_raw[UNCRT==1];

}

transformed parameters {
  vector[N_cases] T_onset[UNCRT];
  vector[n_tot] F[F_is_sampled, sum_D];
  if(UNCRT){
    T_onset[1] = L_ons[1] + (U_ons[1] - L_ons[1]) .* T_raw[1];
  }
  if(F_is_sampled){
#include /chunks/additive_components.stan
  }
}

model {
  
  //Priors
#include /chunks/priors.stan

  // Likelihood model
  if(LH==0){
    
  }else{
    
    if(F_is_sampled){
      
      //  F IS SAMPLED
      // Compute f 
      vector[n_tot] F_sum = rep_vector(0, n_tot);
      for (i in 1:n_tot){
        F_sum[i] += sum(F[1,,i]); 
      }
    
      // Compute likelihood
      if(LH==1){
        // 1. Gaussian likelihood
        y ~ normal(F_sum[1:n], rep_vector(sigma_n[1], n));
      }else{
        real log_g[n] = to_array_1d(F_sum[1:n] + C_hat);
        if(LH==2){
          // 2. Poisson likelihood
          y_int ~ poisson_log(log_g);
        }else if(LH==3){
          // 3. Negative binomial likelihood
          y_int ~ neg_binomial_2_log(log_g, rep_vector(phi[1], n) );
        }else{
          reject("Unknown likelihood!")
        }
      }
    
    }else{
    
      // F NOT SAMPLED
#include /chunks/kernel_matrices.stan
      if(LH!=1){
        reject("Likelihood must be Gaussian if F is not sampled!")
      }
      y  ~ multi_normal_cholesky(mu, Ly);
    }
  }
}

generated quantities {
   vector[n] F_mean_cmp[1 - F_is_sampled, sum_D];
   vector[n] F_var_cmp[ 1 - F_is_sampled, sum_D];
   vector[n] F_mean_tot[1 - F_is_sampled];
   vector[n] F_var_tot[ 1 - F_is_sampled];
                  
   if(F_is_sampled==0){
     matrix[n,n] A;
     vector[n] v;
#include /chunks/kernel_matrices.stan
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

