
<!-- README.md is generated from README.Rmd. Please edit that file -->
hJAM
====

<!-- badges: start -->
<!-- badges: end -->
hJAM is a hierarchical model which unifies the framework of Mendelian Randomization and Transcriptome-wide association studies.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("lailylajiang/hJAM")
```

Example
-------

This is a basic example of fitting hJAM model:

``` r
library(hJAM)
# Download the example data (from simulation) from the package
data(Gl)
data(betas.Gy)
data(conditional_A)
data(marginal_A)
```

If you don't have conditional A matrix, you can use `get_cond_A` (if more than one X) or `get_cond_alpha` (if only one X) to convert the marginal effects to conditional A matrix with the reference panel.

``` r
conditional_A = get_cond_A(marginal_A = marginal_A, Gl = Gl, N.Gx = 1000, ridgeTerm = T)
conditional_alpha = get_cond_alpha(alphas = marginal_A[, 1], Gl = Gl, N.Gx = 1000, ridgeTerm = T)
```

After obtained the conditional A matrix, fit hJAM model with function `hJAM_lnreg`.

``` r
# fit the hJAM model
hJAM_lnreg(betas.Gy = betas.Gy, Gl = Gl, N.Gy = 5000, A = conditional_A, ridgeTerm = T)
#> ------------------------------------------------------ 
#>                    hJAM output                         
#> ------------------------------------------------------ 
#> Number of SNPs used in model: 19 
#> 
#>            Estimate StdErr         95% CI       Pvalue
#> Exposure 1    0.042  0.019 (0.003, 0.082) 3.526250e-02
#> Exposure 2    0.114  0.019 (0.073, 0.155) 2.010226e-05
#> ------------------------------------------------------
```

In the package, you could also implement hJAM with Egger regression, which is designed to detect the unmeasured pleiotropy effect. The function for hJAM with Egger regression is `hJAM_egger`.

``` r
# fit the hJAM model
hJAM_egger(betas.Gy = betas.Gy, Gl = Gl, N.Gy = 5000, A = conditional_A, ridgeTerm = T)
#> ------------------------------------------------------ 
#>                    hJAM egger output                   
#> ------------------------------------------------------ 
#> Number of SNPs used in model: 19 
#> 
#>            Estimate StdErr         95% CI       Pvalue
#> Exposure 1    0.044  0.018 (0.006, 0.082) 2.634095e-02
#> Exposure 2    0.106  0.020 (0.064, 0.148) 6.058969e-05
#> 
#> Intercept
#>      Est.Int  StdErr.Int 95% CI.Int       Pvalue.Int
#> [1,] "-0.895" "0.621"    "(-2.211, 0.42)" "0.168"   
#> ------------------------------------------------------
```
