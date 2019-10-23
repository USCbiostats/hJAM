
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
data(reference_data)
data(betas.Gy)
data(conditional_Z)
data(marginal_Z)
```

If you don't have conditional Z matrix, you can use `get_cond_Z` (if more than one X) or `get_cond_alpha` (if only one X) to convert the marginal effects to conditional Z matrix with the reference panel.

``` r
conditional_Z = get_cond_Z(marginal_Z = marginal_Z, Gl = Gl, N.Gx = 1000, ridgeTerm = T)
conditional_alpha = get_cond_alpha(alphas = marginal_Z[, 1], Gl = Gl, N.Gx = 1000, ridgeTerm = T)
```

After obtained the conditional Z matrix, fit hJAM model with function `hJAM_lnreg`.

``` r
# fit the hJAM model
hJAM_lnreg(betas.Gy = betas.Gy, Gl = Gl, N.Gy = 5000, Z = conditional_Z, ridgeTerm = T)
#> -------------------------------------------- 
#>                 hJAM output                  
#> -------------------------------------------- 
#> Number of SNPs used in model: 19 
#> 
#>              Estimate     StdErr       Pvalue
#> Exposure 1 0.04243959 0.01855368 3.526250e-02
#> Exposure 2 0.11365449 0.01949630 2.010226e-05
#> --------------------------------------------
```

In the package, you could also implement hJAM with Egger regression, which is designed to detect the unmeasured pleiotropy effect. The function for hJAM with Egger regression is `hJAM_egger`.

``` r
# fit the hJAM model
hJAM_egger(betas.Gy = betas.Gy, Gl = Gl, N.Gy = 5000, Z = conditional_Z, ridgeTerm = T)
#> -------------------------------------------- 
#>              hJAM egger output               
#> -------------------------------------------- 
#> Number of SNPs used in model: 19 
#> 
#> Exposures
#>              Estimate     StdErr       Pvalue
#> Exposure 1 0.04411007 0.01802767 2.634095e-02
#> Exposure 2 0.10587960 0.01965768 6.058969e-05
#> 
#> Intercept
#>         Est.Int StdErr.Int Pvalue.Int
#> [1,] -0.8951416  0.6205016  0.1684211
#> --------------------------------------------
```
