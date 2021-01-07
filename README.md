
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hJAM

<!-- badges: start -->

[![R build
status](https://github.com/lailylajiang/hJAM/workflows/R-CMD-check/badge.svg)](https://github.com/lailylajiang/hJAM)
<!-- badges: end -->

<!-- CRAN badges: start -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/hJAM)](https://CRAN.R-project.org/package=hJAM)[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/hJAM)](https://CRAN.R-project.org/package=hJAM)
<!-- CRAN badges: end -->

hJAM is a hierarchical model which unifies the framework of Mendelian
Randomization and Transcriptome-wide association studies.

## Installation

You can install the published version of hJAM from CRAN with:

``` r
install.packages("hJAM")
```

Or you can install the development version from
[GitHub](https://github.com/lailylajiang/hJAM) with:

``` r
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("lailylajiang/hJAM")
```

## Example

This is a basic example of fitting hJAM model:

``` r
library(hJAM)
#> 
# Download the data for data example 2 from the package
data(MI)
```

## Data and functions import

A quick look at the data in the example -

``` r
MI.Amatrix[1:5, ]
#>               bmi          t2d
#> [1,]  0.019531085  0.072587211
#> [2,]  0.025262061  0.013586392
#> [3,] -0.005147363  0.089673178
#> [4,]  0.046302578  0.041313103
#> [5,]  0.016849395 -0.004564683
MI.betas.gwas[1:5]
#> [1]  0.0197298348  0.0133151413  0.0008717583  0.0213550014 -0.0031514278
MI.SNPs_info[1:5, ]
#>          SNP Major_A   ref_frq BMI.sig T2D.sig
#> 1  rs2296173       G 0.1207945       0       1
#> 2   rs657452       A 0.5474260       1       0
#> 3 rs12088739       A 0.8785975       0       1
#> 4  rs3101336       C 0.6781516       1       0
#> 5 rs12566985       G 0.6793677       1       0
```

``` r
MI.cond_A = JAM_A(marginalA = MI.marginal.Amatrix, Geno = MI.Geno, N.Gx = c(339224, 659316), ridgeTerm = TRUE)
MI.cond_A[1:5, ]
#>               bmi          t2d
#> [1,]  0.019531085  0.072587479
#> [2,]  0.025262061  0.013586305
#> [3,] -0.005147363  0.089673869
#> [4,]  0.046302578  0.041313424
#> [5,]  0.016849395 -0.004564612
MI.Amatrix[1:5, ]
#>               bmi          t2d
#> [1,]  0.019531085  0.072587211
#> [2,]  0.025262061  0.013586392
#> [3,] -0.005147363  0.089673178
#> [4,]  0.046302578  0.041313103
#> [5,]  0.016849395 -0.004564683
```

``` r
hJAM::hJAM(betas.Gy = MI.betas.gwas, N.Gy = 459324, A = MI.Amatrix, 
           Geno = MI.Geno, ridgeTerm = TRUE) # 459324 is the sample size of the UK Biobank GWAS of MI
#> ------------------------------------------------------ 
#>                    hJAM output                         
#> ------------------------------------------------------ 
#> Number of SNPs used in model: 210 
#> 
#>     Estimate StdErr         95% CI       Pvalue
#> bmi    0.322  0.061 (0.202, 0.441) 1.268210e-07
#> t2d    0.119  0.017 (0.086, 0.153) 3.176604e-12
#> ------------------------------------------------------
```
