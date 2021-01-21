
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hJAM

<!-- badges: start -->

[![R build
status](https://github.com/lailylajiang/hJAM/workflows/R-CMD-check/badge.svg)](https://github.com/lailylajiang/hJAM)
<!-- badges: end -->

<!-- CRAN badges: start -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/hJAM)](https://CRAN.R-project.org/package=hJAM)[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/hJAM)](https://CRAN.R-project.org/package=hJAM)
<!-- CRAN badges: end -->

## About

In this package, we provide three methods that are developed for
estimating (and selecting) the causal effect of intermediates on the
outcome using the genetic variants as instrumental variables.

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

## Citing this work

If you find the `hJAM`, `hJAM_egger`, and/or `JAM_A` useful, please
cite:

  - Jiang, L., Xu, S., Mancuso, N., Newcombe, P.J. & Conti, D.V. A
    Hierarchical Approach Using Marginal Summary Statistics for Multiple
    Intermediates in a Mendelian Randomization or Transcriptome
    Analysis. (Accepted by American Journal of Epidemiology).

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
#> R version 4.0.2 (2020-06-22)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Catalina 10.15.7
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.0.2  magrittr_2.0.1  tools_4.0.2     htmltools_0.5.0
#>  [5] yaml_2.2.1      stringi_1.5.3   rmarkdown_2.3   knitr_1.30     
#>  [9] stringr_1.4.0   xfun_0.20       digest_0.6.27   rlang_0.4.10   
#> [13] evaluate_0.14
```
