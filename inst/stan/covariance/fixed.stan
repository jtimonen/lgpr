{
  int ix;
  KF[1] = STAN_K_zerosum(X[1], X[1], N_cat[1]);
  for(j in 1:D[3]){
    KF[1+j] = STAN_K_bin(X_notnan, X_notnan, 1);
  }
  for(j in 1:D[5]){
    ix = 2 + D[3] + D[4] + j;
    KF[1+D[3]+j] = STAN_K_zerosum(X[ix], X[ix], N_cat[1+j]);
  }
  for(j in 1:D[6]){
    ix = 2 + D[3] + D[4] + D[5] + j;
    KF[1+D[3]+D[5]+j] = STAN_K_zerosum(X[ix], X[ix], N_cat[1+D[5]+j]);
  }
}
