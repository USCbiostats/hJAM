
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

After obtained the conditional Z matrix, fit hJAM model with function `hJAM`.

``` r
# fit the hJAM model
hJAM(betas.Gy = betas.Gy, Gl = Gl, N.Gy = 5000, Z = conditional_Z, ridgeTerm = T)
#> $Estimate
#>         X1         X2 
#> 0.06600536 0.03975728 
#> 
#> $StdErr
#>          X1          X2 
#> 0.009699963 0.008280250 
#> 
#> $Pvalue
#>           X1           X2 
#> 2.263574e-06 1.429418e-04
```
