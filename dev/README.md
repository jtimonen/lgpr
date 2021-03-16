# lgpr/dev

These files are only for development purposes and not a part of the package
itself. They are supposed to be run while in the package root directory.

* `cpp.R` builds the Stan functions into C++ code
* `covr.R` creates a test coverage report
* `exports.R` lists the exported functions
* `style.R` corrects styling errors the package code
* `lint.R` finds styling errors in the package code


## NOTES (1.1.9031)

```
installing to ...
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
  converting help for package 'lgpr'
    finding HTML links ... done
  
    colorset                                html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/colorset.Rd:11: file link 'color_scheme_get' in package 'bayesplot' does not exist and so has been treated as a topic

    lgp                                     html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/lgp.Rd:89: file link 'sampling' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/lgp.Rd:89: file link 'optimizing' in package 'rstan' does not exist and so has been treated as a topic

    lgpfit-class                            html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/lgpfit-class.Rd:53: file link 'stanfit' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/lgpfit-class.Rd:57: file link 'stanfit' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/lgpfit-class.Rd:59: file link 'extract' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/lgpfit-class.Rd:63: file link 'extract' in package 'rstan' does not exist and so has been treated as a topic
 
    plot_draws                              html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/plot_draws.Rd:38: file link 'mcmc_intervals' in package 'bayesplot' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/plot_draws.Rd:38: file link 'mcmc_dens' in package 'bayesplot' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/plot_draws.Rd:39: file link 'mcmc_areas' in package 'bayesplot' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/plot_draws.Rd:39: file link 'mcmc_trace' in package 'bayesplot' does not exist and so has been treated as a topic


    sample_model                            html  
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/sample_model.Rd:21: file link 'sampling' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/sample_model.Rd:24: file link 'optimizing' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/sample_model.Rd:16: file link 'sampling' in package 'rstan' does not exist and so has been treated as a topic
Rd warning: C:/Users/Juho/AppData/Local/Temp/RtmpUJ8vRJ/R.INSTALL2e8013ea22c9/lgpr/man/sample_model.Rd:16: file link 'optimizing' in package 'rstan' does not exist and so has been treated as a topic


** building package indices
** testing if installed package can be loaded from temporary location
*** arch - i386
*** arch - x64
** testing if installed package can be loaded from final location
*** arch - i386
*** arch - x64
** testing if installed package keeps a record of temporary installation path
* DONE (lgpr)
> library(lgpr)
Attached lgpr 1.1.0.9031, using rstan 2.21.2.

> ?lgpr
> a <- simulate_data(N = 8, seq(0,100,by=2))
> a
An object of class lgpsim. See ?lgpsim for more info.
> plot(a)
Error in as.double(y) : 
  cannot coerce type 'S4' to vector of type 'double'
```

