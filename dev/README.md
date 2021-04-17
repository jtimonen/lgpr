# lgpr/dev

These files are only for development purposes and not a part of the package
itself. They are supposed to be run while in the package root directory.

* `covr.R` creates a test coverage report
* `exports.R` lists the exported functions
* `style.R` corrects styling errors the package code
* `lint.R` finds styling errors in the package code
* (`cpp.R` builds the Stan functions into C++ code)

## NOTES (1.1.9000)

```
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/kernels_fpost.Rd:39: file link 'get_stream' in package 'rstan' does not exist and so has been treated as a topic
    lgp                                     html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/lgp.Rd:93: file link 'sampling' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/lgp.Rd:100: file link 'sampling' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/lgp.Rd:100: file link 'optimizing' in package 'rstan' does not exist and so has been treated as a topic
    lgpexpr-class                           html  
    lgpfit-class                            html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/lgpfit-class.Rd:80: file link 'stanfit' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/lgpfit-class.Rd:84: file link 'stanfit' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/lgpfit-class.Rd:89: file link 'extract' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/lgpfit-class.Rd:93: file link 'extract' in package 'rstan' does not exist and so has been treated as a topic
    lgpformula-class                        html  
    lgpmodel-class                          html  
    lgpr-package                            html  
    lgprhs-class                            html  
    lgpscaling-class                        html  
    lgpsim-class                            html  
    lgpterm-class                           html  
    model_summary                           html  
    new_x                                   html  
    operations                              html  
    plot_api_c                              html  
    plot_api_g                              html  
    plot_components                         html  
    plot_data                               html  
    plot_draws                              html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/plot_draws.Rd:38: file link 'mcmc_intervals' in package 'bayesplot' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/plot_draws.Rd:38: file link 'mcmc_dens' in package 'bayesplot' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/plot_draws.Rd:39: file link 'mcmc_areas' in package 'bayesplot' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/plot_draws.Rd:39: file link 'mcmc_trace' in package 'bayesplot' does not exist and so has been treated as a topic
    plot_inputwarp                          html  
    plot_invgamma                           html  
    plot_pred                               html  
    plot_sim                                html  
    posterior_f                             html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/posterior_f.Rd:38: file link 'get_stream' in package 'rstan' does not exist and so has been treated as a topic
    ppc                                     html  
    pred                                    html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/pred.Rd:37: file link 'get_stream' in package 'rstan' does not exist and so has been treated as a topic
    prior_to_num                            html  
    priors                                  html  
    relevances                              html  
    sample_model                            html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/sample_model.Rd:37: file link 'sampling' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/sample_model.Rd:40: file link 'optimizing' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/sample_model.Rd:25: file link 'sampling' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/sample_model.Rd:32: file link 'sampling' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpuA0e7t/R.INSTALL1560350a78e/lgpr/man/sample_model.Rd:32: file link 'optimizing' in package 'rstan' does not exist and so has been treated as a topic
    select                                  html  
    show                                    html  
    sim.create_f                            html  
    sim.create_x                            html  
    sim.create_y                            html  
    sim.kernels                             html  
    simulate_data                           html  
    split                                   html  
    testdata_001                            html  
    testdata_002                            html  
    validate                                html  
    var_mask                                html  
    warp_input                              html  
```

