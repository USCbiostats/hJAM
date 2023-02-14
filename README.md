
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JAM: Joint Analysis of Marginal Summary Statistics

<!-- badges: start -->

[![R build
status](https://github.com/lailylajiang/hJAM/workflows/R-CMD-check/badge.svg)](https://github.com/lailylajiang/hJAM)
<!-- badges: end -->

<!-- CRAN badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/hJAM)](https://CRAN.R-project.org/package=hJAM)[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/hJAM)](https://CRAN.R-project.org/package=hJAM)
<!-- CRAN badges: end -->

## About

In this package, we provide three methods that utilizes GWAS summary
statistics to perform post-GWAS analysis, including Mendelian
Randomization, TWAS, and multi-population fine-mapping.

## hJAM and SHAJAM

The hJAM (`hJAM`) is a hierarchical model which unifies the framework of
Mendelian Randomization and Transcriptome-wide association studies. The
hJAM-Egger (`hJAM_egger`) is a natural extension on hJAM that uses an
intercept term to account for the pleiotropy effect of the SNPs on the
outcome. The SHA-JAM (`SHAJAM`) is applicable for high-throughput
experiment data, such as omics data.

Additionally, we provide implementations to construct the weight matrix.
The JAM framework (`JAM_A`) is used to convert the marginal summary
statistics into conditional ones by using the correlation matrix of the
SNPs from a reference data. The SuSiE JAM (`susieJAM_A`) is used to
select the SNPs for one intermediate using the marginal summary
statistics from GWAS or other study summary data. It can also be used
for fine-mapping problems.

## mJAM

mJAM is for multi-population fine-mapping using GWAS summary statistics
and credible set construction. A tutorial to get started with mJAM can
be found [here](https://jiayi-s.github.io/more_on_mJAM/).

## Citing this work

If you find the `hJAM`, `hJAM_egger`, and/or `JAM_A` useful, please
cite:

- Jiang, L., Xu, S., Mancuso, N., Newcombe, P.J. & Conti, D.V. A
  Hierarchical Approach Using Marginal Summary Statistics for Multiple
  Intermediates in a Mendelian Randomization or Transcriptome Analysis.
  (Accepted by American Journal of Epidemiology).

If you find the `SHAJAM` and/or `susieJAM_A` useful, please cite:

- Jiang, L., Conti, D.V. SHA-JAM: A Scalable Hierarchical Approach to
  Joint Analysis for Marginal Summary Statistics with Omics Data. (In
  preparation; preprint will be ready in Feburary, 2021).

## Quick start with hJAM package

You can install the published version of hJAM from CRAN with:

``` r
install.packages("hJAM")
```

Or you can install the development version from
[GitHub](https://github.com/USCbiostats/hJAM) with:

``` r
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("USCbiostats/hJAM")
```

## Session info

``` r
sessionInfo()
#> R version 4.2.0 (2022-04-22)
#> Platform: aarch64-apple-darwin20 (64-bit)
#> Running under: macOS 13.0.1
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.2.0  magrittr_2.0.3  fastmap_1.1.0   cli_3.4.0      
#>  [5] tools_4.2.0     htmltools_0.5.3 rstudioapi_0.14 yaml_2.3.5     
#>  [9] stringi_1.7.8   rmarkdown_2.14  knitr_1.40      stringr_1.4.1  
#> [13] xfun_0.32       digest_0.6.29   rlang_1.0.5     evaluate_0.16
```
