vector[N_cases] T_effect[UNCRT];
if(UNCRT){
  T_effect[1] = L_ons[1] + (U_ons[1] - L_ons[1]) .* T_raw[1];
}
