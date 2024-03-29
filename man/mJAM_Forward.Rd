% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mJAM_Forward.R
\name{mJAM_Forward}
\alias{mJAM_Forward}
\title{Run mJAM with Forward Selection}
\usage{
mJAM_Forward(
  N_GWAS,
  X_ref,
  Marg_Result,
  EAF_Result,
  condp_cut = NULL,
  index_snps = NULL,
  within_pop_threshold = 0.5,
  across_pop_threshold = 0.2,
  coverage = 0.95,
  Pr_Med_cut = 0,
  filter_rare = FALSE,
  rare_freq = NULL,
  filter_unstable_est = FALSE,
  use_robust_var_est = FALSE
)
}
\arguments{
\item{N_GWAS}{A vector of sample sizes in all original GWAS studies.}

\item{X_ref}{A list of matrices with individual-level SNP dosage data in each study/population. Each column corresponds to a SNP. Note that the columns name should match exactly to the SNP column in `Marg_Result` and `EAF_Result`. If certain SNP(s) is missing in dosage, then insert NAs in corresponding column(s).}

\item{Marg_Result}{A data frame with marginal summary statistics from all studies. Col1: SNP name; Col2: Effect sizes from study #1; Col3: Std Errors of effect sizes from study #1; ...}

\item{EAF_Result}{A data frame with effect allele frequency (EAF) from all studies. Col1: SNP name; Col2: EAF from study #1; Col3: EAF from study #2; ...}

\item{condp_cut}{Threshold of conditional p-value to be considered as significant. No default specified. Usually recommend 5e-8.}

\item{index_snps}{User-defined index SNP(s), if any. Default is `NULL` which means mJAM-Forward will automatically select index variants.}

\item{within_pop_threshold}{Threshold of r2 with selected index SNP(s) within a single population. If a SNP's correlation with any selected index SNP is greater than this threshold in at least one population, it will be excluded from subsequent rounds of index SNP selection.}

\item{across_pop_threshold}{Threshold of r2 with selected index SNP(s) across all populations. If a SNP's correlation with any selected index SNP is greater than this threshold in all populations, it will be excluded from subsequent rounds of index SNP selection.}

\item{coverage}{The required coverage of credible sets. Default is 0.95.}

\item{Pr_Med_cut}{Cut off of mJAM posterior mediation probability (P(Med)) during credible set construction. Low P(Med) may indicate low correlation between the candidate SNP and the index SNP. Any candidate credible set SNPs with P(Med) < Pr_Med_cut will be not be considered for credible set. Default is 0.}

\item{filter_rare}{A logical variable indicating whether to filter rare SNPs before the analysis. Default is `FALSE.` If `TRUE`, then please specify `rare_freq`.}

\item{rare_freq}{A vector of frequencies between 0 and 0.5 to specify the minor allele frequency cut-off if you want to filter rare SNPs before the analysis. Please also set `filter_rare` to be TRUE. For example, if there are 3 populations, then rare_freq = c(0.01, 0, 0.01) means SNPs with MAF < 0.01 in pop 1 and MAF < 0.01 in pop 3 will be removed from analysis.}

\item{filter_unstable_est}{whether to filter variants with inconsistent estimate between mJAM and meta-analysis.}

\item{use_robust_var_est}{whether to use the robust estimate of residual variance (weighting between median and individual estimates).}
}
\value{
\describe{
   \item{index}{A table listing all the selected index SNP(s) (`SNP`), along with their log10(p-value) conditional on all SNP(s) above (`cond_log10p`), the log10(p-value) conditional on all other index SNP(s) (`final_log10p`), and the p-value threshold used in this analysis (`pcut`).}
   \item{cs}{A table recording various posterior probabilities of all SNPs being considered for credible set SNPs. }
   \item{mJAM_marg_est}{A table with the marginal effect estimates and standard errors of all SNPs under the mJAM model.}
   \item{QC_marg_est}{The complete table of marginal effect estimates using fixed-effect model and mJAM model. For QC purpose only.}
}
}
\description{
fitting mJAM-Forward
}
\author{
Jiayi Shen
}
