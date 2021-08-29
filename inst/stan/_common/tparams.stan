  // Transform x uncertainty parameters
  vector[num_xpar] xpar[num_unc>0];
  if(num_unc > 0) {
    xpar[1] = xpar_lb[1] + (xpar_ub[1] - xpar_lb[1]) .* xpar_raw[1];
  }
