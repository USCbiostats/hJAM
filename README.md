
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
Randomization and Transcriptome-wide association studies. This package
contains the implementations of the methods that were proposed in two
papers from our lab:

  - Jiang, L., Xu, S., Mancuso, N., Newcombe, P.J. & Conti, D.V. A
    Hierarchical Approach Using Marginal Summary Statistics for Multiple
    Intermediates in a Mendelian Randomization or Transcriptome
    Analysis. (Accepted by American Journal of Epidemiology).
  - Jiang, L., Conti, D.V. SHA-JAM: A Scalable Hierarchical Approach to
    Joint Analysis for Marginal Summary Statistics with Omics Data. (In
    preparation).

Please cite if you find the `hJAM` package or any of the source code in
this repository useful for your study.

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
