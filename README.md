
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JAM: Joint Analysis of Marginal summary statistics

<!-- badges: start -->

[![R build
status](https://github.com/lailylajiang/hJAM/workflows/R-CMD-check/badge.svg)](https://github.com/lailylajiang/hJAM)
<!-- badges: end -->

<!-- CRAN badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/hJAM)](https://CRAN.R-project.org/package=hJAM)[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/hJAM)](https://CRAN.R-project.org/package=hJAM)
<!-- CRAN badges: end -->

## About

In this package, we provide three `JAM`-family methods that utilize GWAS
summary statistics to perform post-GWAS analysis, including Mendelian
Randomization, TWAS, and multi-population fine-mapping.

## hJAM

hJAM (`hJAM`) is a hierarchical model which unifies the framework of
Mendelian Randomization and Transcriptome-wide association studies. The
hJAM-Egger (`hJAM_egger`) is a natural extension on hJAM that uses an
intercept term to account for the pleiotropy effect of the SNPs on the
outcome.

Additionally, we provide implementations to construct the weight matrix
`A` in hJAM. `JAM_A` converts the marginal summary statistics into
conditional ones by using the correlation matrix of the SNPs from a
reference data. `susieJAM_A` selects SNPs for one intermediate using the
marginal summary statistics from GWAS or other study summary data.

## SHA-JAM

SHA-JAM (`SHAJAM`) is a scalable version of `hJAM` that handles
high-dimensional intermediates. SHA-JAM performs model selection from
highly correlated intermediates through SuSiE: Sum of Single Effect
Model (`SHAJAM`) or elastic-net (`EN.hJAM`).

## mJAM

mJAM is for multi-population fine-mapping using GWAS summary statistics
and credible set construction. We provide two implementations of mJAM:
one through SuSiE (`mJAM_SuSiE`) and another one through forward
selection (`mJAM_Forward`). `mJAM_Forward` also provides the flexiblity
of constructing credible sets for using user-defined index SNPs. A
tutorial to get started with mJAM can be found
[here](https://jiayi-s.github.io/more_on_mJAM/).

## Citing this work

- **hJAM**: Jiang, L., Xu, S., Mancuso, N., Newcombe, P. J., &
  Conti, D. V. (2021). A Hierarchical Approach Using Marginal Summary
  Statistics for Multiple Intermediates in a Mendelian Randomization or
  Transcriptome Analysis. *American journal of epidemiology*, 190(6),
  1148–1158. <https://doi.org/10.1093/aje/kwaa287>

- **SHA-JAM**: Jiang, L., Conti, D.V. SHA-JAM: A Scalable Hierarchical
  Approach to Joint Analysis for Marginal Summary Statistics with Omics
  Data. (In preparation).

- **mJAM**: Shen, J., Jiang, L., Wang, K., Wang, A., Chen, F., Newcombe,
  P.J., Haiman, C.A., & Conti, D.V. Fine-Mapping and Credible Set
  Construction using a Multi-population Joint Analysis of Marginal
  Summary Statistics from Genome-wide Association Studies.
  ([Preprint](https://www.biorxiv.org/content/10.1101/2022.12.22.521659v1)
  available)

## Quick start with hJAM package

You can install the published version of hJAM from CRAN with:

``` r
install.packages("hJAM")
```

Currently we are working on improving `hJAM` package and adding
genome-wide implementation of mJAM. If you want to use the most updated
version of this package, we recommend installing the development version
from [GitHub](https://github.com/USCbiostats/hJAM) with:

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
