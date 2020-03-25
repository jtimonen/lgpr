// D is an array of six integers, so that
//   D[1] = binary value indicating if f(id,age) is in model
//   D[2] = binary value indicating if f(age) is in model
//   D[3] = binary value indicating if f(diseaseAge) is in model
//   D[4] = #{other continuous covariates}
//   D[5] = #{discr. covars. that cause deviation from f(age)}
//   D[6] = #{discr. covars. that only have a baseline offset effect}
int<lower=0> D[6];

int<lower=1> N_tot;     // #{individuals}
int<lower=0> N_cases;   // #{case individuals}
int<lower=2> d;         // #{covariates} (id and age required)
int<lower=1> n;         // #{observations}
vector[n] X[d];         // covariates, X[j] is the jth covariate
int       X_notnan[n];  // X_notnan[i] tells if X_diseaseAge[i] isn't NaN
vector[n] y;            // the response variable (as a vector of reals)

// #{categs} for each categorical covariate
int<lower=1> N_cat[1+D[5]+D[6]]; 

// Option switches
int<lower=0,upper=1> VERBOSE;      // is model info printed?
int<lower=0,upper=1> USE_VAR_MASK; // is variance mask kernel used?
int<lower=0,upper=1> UNCRT;        // are diseaseAge measurements uncertain?
int<lower=0,upper=1> HMGNS;        // is diseaseAge effect is homogenous?
int<lower=0,upper=1> BACKWARDS;    // is prior of effect time "backwards"?
int<lower=0,upper=1> RELATIVE;     // is prior of effect time rel. to obs. one?
int<lower=0,upper=1> SKIP_LH;      // skip likelihood evaluation?
int<lower=0,upper=1> SKIP_GQ;      // skip the generated quantities block?

// Prior types and transforms
int t_ID[D[1],4];         // for id*age component
int t_A[D[2],4];          // for shared age component
int t_D[D[3],6];          // for disease age component
int t_CNT[D[4],4];        // for components with continuous covariates
int t_CAT[D[5],4];        // for components with discrete covariates
int t_OFS[D[6],2];        // for offset components
int t_ONS[N_cases,2];     // for onset, if uncertain
  
// Hyperparameters of the priors and scaling factors,
real p_ID[D[1],6];         // for id*age component
real p_A[D[2],6];          // for shared age component
real p_D[D[3],9];          // for disease age component
real p_CNT[D[4],6];        // for components with continuous covariates
real p_CAT[D[5],6];        // for components with discrete covariates
real p_OFS[D[6],3];        // for offset components
real p_BET[2];             // hyperparameters of the beta prior of BETA
real p_ONS[N_cases, 3];    // for onset, if uncertain
  
// Inputs related to mapping from row index to case index and back
int<lower=0,upper=n>        M_max;
int<lower=0>                caseID_to_rows[N_cases, M_max];
int<lower=0,upper=N_cases>  row_to_caseID_plus1[n];
int<lower=0,upper=M_max>    caseID_nrows[N_cases];

// Inputs related to uncertain disease onset
vector[N_cases] T_observed;   // observed disease effect times
vector[N_cases] L_ons[UNCRT]; // low bounds for disease effect times
vector[N_cases] U_ons[UNCRT]; // upper bounds for disease effect times
  
// Other
real DELTA;         // jitter to ensure pos. def. kernel matrices
real vm_params[2];  // variance mask parameters
