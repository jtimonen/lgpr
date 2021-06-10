# lgpr/dev

These files are only for development purposes and not a part of the package
itself. They are supposed to be run while in the package root directory.

* `covr.R` creates a test coverage report
* `exports.R` lists the exported functions
* `style.R` corrects styling errors the package code
* `lint.R` finds styling errors in the package code
* (`cpp.R` builds the Stan functions into C++ code)

## Test suite output

* Version: `1.1.9000`
* Commit: `1a84c1bac5e5263e259f4b2bb58ac7978c8a205`
* Date: `9th June 2021`
* TotalCPUTime: `04:54:31`
* MaxRSS: `3240M`

Proteomics data experiment's (4-5) relevance results confirmed to agree with paper.

```

                   name obs_model f_sampled n_obs n_comps t_ch_mean t_ch_sd
1         01_sleepstudy  gaussian     FALSE   180       2    170.74    3.01
2          02_orthodont  gaussian     FALSE   108       3     74.94    4.77
3 03_two-age-components  gaussian     FALSE    70       2     14.07    0.85
4        04_protein-liu  gaussian     FALSE   189       5    524.67   20.78
5  05_protein-liu-heter  gaussian     FALSE   189       5   1360.18   59.49

     t_fit t_pred t_post  t_total size_fit size_small size_disk n_div max_Rhat
1  1397.94  44.95   9.58  1456.06  44.9 Mb       1 Mb   40.1 Mb     0    1.002
2   615.33  23.24  10.00   650.88  34.1 Mb     1.1 Mb   27.9 Mb     0    1.002
3   119.50  11.55   3.83   136.21  18.1 Mb       1 Mb   17.3 Mb     0    1.004
4  4270.24 102.16  15.64  4393.96  82.1 Mb     1.4 Mb   63.6 Mb     0    1.005
5 10966.25  97.41  17.81 11086.60  82.8 Mb     2.1 Mb   64.5 Mb     0    1.005

  min_TESS min_BESS
1     2517     2036
2     2087     1753
3     2066     1662
4     1344     1470
5     1015      690


-----------------------------------------------------------------------
01_sleepstudy.R: 
           gp(Days) gp(Days)*zs(Subject)    noise
Relevance 0.3096018            0.5637252 0.126673
Expected  0.3096018            0.5637252 0.126673
-----------------------------------------------------------------------
02_orthodont.R: 
            gp(age) gp(age)*zs(Subject) gp(age)*zs(Sex)     noise
Relevance 0.2818898           0.4012653       0.1506231 0.1662218
Expected  0.2818898           0.4012653       0.1506231 0.1662218
-----------------------------------------------------------------------
03_two-age-components.R: 
            gp(age) gp(age_fast)     noise
Relevance 0.1034861    0.7895722 0.1069417
Expected  0.0764848    0.8130515 0.1104637
-----------------------------------------------------------------------
04_protein-liu.R: 
          gp(age)*zs(id)   gp(age) gp_vm(diseaseAge) zs(sex)*gp(age)
Relevance      0.2286768 0.1573663        0.02842119      0.03088412
Expected       0.2286768 0.1573663        0.02842119      0.03088412
            zs(group)     noise
Relevance 0.007530387 0.5471212
Expected  0.007530387 0.5471212
-----------------------------------------------------------------------
05_protein-liu-heter.R: 
          gp(age)*zs(id)   gp(age) het(id)*gp_vm(diseaseAge) zs(sex)*gp(age)
Relevance     0.09591313 0.1149214                 0.2531681      0.03665491
Expected      0.09591313 0.1149214                 0.2531681      0.03665491
            zs(group)     noise
Relevance 0.003777137 0.4955653
Expected  0.003777137 0.4955653
-----------------------------------------------------------------------


R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)
```

