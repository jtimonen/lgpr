  // Transform raw effect times
  vector[num_teff] teff[idx_unc>0];
  if(idx_unc > 0) {
    teff[1] = teff_lb[1] + (teff_ub[1] - teff_lb[1]) .* teff_raw[1];
  }
