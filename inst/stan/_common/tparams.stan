  // Transform raw effect times
  array[num_uncrt>0] vector[num_bt] teff;
  for(j in 1:num_uncrt){
    teff[j] = teff_lb[j] + (teff_ub[j] - teff_lb[j]) .* teff_raw[j];
  }
