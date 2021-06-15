# lgpr/dev

These files are only for development purposes and not a part of the package
itself. They are supposed to be run while in the package root directory.

* `covr.R` creates a test coverage report
* `exports.R` lists the exported functions
* `style.R` corrects styling errors the package code
* `lint.R` finds styling errors in the package code
* (`cpp.R` builds the Stan functions into C++ code)

## Test suite output

* Version: `1.1.0`
* Date: `13th June 2021`
* TotalCPUTime: `05:27:25`
* MaxRSS: `3524M`

Proteomics data experiment's (4-5) relevance results confirmed to agree with paper.

```
                   name obs_model f_sampled n_obs n_comps t_ch_mean t_ch_sd
1         01_sleepstudy  gaussian     FALSE   180       2    193.56    3.13
2          02_orthodont  gaussian     FALSE   108       3     78.50    5.00
3 03_two-age-components  gaussian     FALSE    70       2     16.10    0.23
4        04_protein-liu  gaussian     FALSE   189       5    572.94   22.37
5  05_protein-liu-heter  gaussian     FALSE   189       5   1528.20   66.63

     t_fit t_pred t_post  t_total size_fit size_small size_disk n_div max_Rhat
1  1578.72  51.27  11.23  1645.33  44.9 Mb       1 Mb   40.1 Mb     0    1.002
2   643.12  24.99  10.82   681.67  34.1 Mb     1.1 Mb   27.9 Mb     0    1.002
3   136.60  12.39   4.16   154.59  18.1 Mb       1 Mb   17.3 Mb     0    1.002
4  4642.29  98.22  16.55  4763.23  82.1 Mb     1.4 Mb   63.6 Mb     0    1.005
5 12283.78 100.55  18.40 12408.22  82.8 Mb     2.1 Mb   64.5 Mb     0    1.005

  min_TESS min_BESS
1     2517     2036
2     2087     1753
3     1775     1618
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
Relevance 0.05237128    0.8587906 0.0888381
Expected  0.07648480    0.8130515 0.1104637

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


                     file                                          rds
1         01_sleepstudy.R         test_suite_out/rds/01_sleepstudy.rds
2          02_orthodont.R          test_suite_out/rds/02_orthodont.rds
3 03_two-age-components.R test_suite_out/rds/03_two-age-components.rds
4        04_protein-liu.R        test_suite_out/rds/04_protein-liu.rds
5  05_protein-liu-heter.R  test_suite_out/rds/05_protein-liu-heter.rds


R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)
```

