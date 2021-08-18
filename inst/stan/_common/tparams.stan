  // Transform raw effect times
  vector[num_teff] teff[num_unc>0];
  if(num_unc > 0) {
    teff[1] = teff_lb[1] + (teff_ub[1] - teff_lb[1]) .* teff_raw[1];
  }
