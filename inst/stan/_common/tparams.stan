  // Transform raw effect times
  vector[num_bt] teff[num_uncrt>0];
  for(j in 1:num_uncrt){
    teff[j] = teff_lb[j] + (teff_ub[j] - teff_lb[j]) .* teff_raw[j];
  }
