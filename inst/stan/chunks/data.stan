  
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

// Kernel types that have options
int<lower=0,upper=1> USE_VAR_MASK;
int<lower=0,upper=1> cat_interact_kernel; // {1 = categorical, 0 = binary} kernel
real vm_params[2];

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
                  
// Verbosity level
int<lower=0,upper=1> verbose_mode;

