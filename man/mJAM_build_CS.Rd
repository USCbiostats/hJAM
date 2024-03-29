% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mJAM_build_CS.R
\name{mJAM_build_CS}
\alias{mJAM_build_CS}
\title{Construct mJAM credible set based for selected index SNP}
\usage{
mJAM_build_CS(
  X_id,
  prev_X_list = NULL,
  All_id,
  PrCS_weights = "Pr(M_C)",
  coverage = 0.95,
  GItGI_curr,
  GIty_curr,
  yty_curr,
  yty_med,
  N_GWAS,
  rare_SNPs = NULL,
  Pr_Med_cut = 0.1,
  use_robust_var_est = FALSE
)
}
\arguments{
\item{X_id}{A character specifying the ID of the index SNP; should be found in `All_id`.}

\item{prev_X_list}{A list of character vector of the ID(s) of previously selected index SNP(s).}

\item{All_id}{A list of character vector of the ID(s) of all SNP(s) remaining in the analysis, including all previously selected SNP(s) and the current index SNP.}

\item{PrCS_weights}{An option to specify what weights to apply on Pr(Med). Default is "Pr(M_C)".}

\item{coverage}{A number between 0 and 1 specifying the “coverage” of the estimated confidence sets.}

\item{GItGI_curr}{A list of GItGI statistics at the current stage (after pruning out SNPs correlated with previously selected index SNPs).}

\item{GIty_curr}{A list of GIty estimates of all remaining SNPs at the current stage (after pruning out SNPs correlated with previously selected index SNPs).}

\item{yty_curr}{A list of yty estimates of all remaining SNPs at the current stage (after pruning out SNPs correlated with previously selected index SNPs).}

\item{yty_med}{A list of median yty across all SNPs.}

\item{N_GWAS}{A vector of sample sizes in all original GWAS studies.}

\item{rare_SNPs}{A numeric vector of ID(s) for rare SNP(s) which we do not apply weighting. Instead, we use the individual estimate of yty for these SNPs for robustness.}

\item{Pr_Med_cut}{The cutoff for Pr(Mediation); SNPs with Pr(Mediation) smaller than this cutoff will be assigned a Pr(CS) = 0 and thus not included in the credible set for the current index}

\item{use_robust_var_est}{whether to use linear combination of median yty and individual yty.}
}
\value{
A table with the following columns:

\describe{
   \item{CS_SNP}{SNP name.}
   \item{Post_Model_Prob}{The posterior Pr(Model) of this SNP on its absolute scale.}
   \item{Post_Model_Prob_Ratio}{The posterior Pr(Model) of this SNP divided by the posterior Pr(Model) of index SNP. It should be <= 1.}
   \item{Post_Model_Prob_Ratio2}{If `Post_Model_Prob_Ratio` is greater than 1, set `Post_Model_Prob_Ratio2` to 1. Otherwise, it's the same as `Post_Model_Prob_Ratio`.}
   \item{Med_Effect_Size}{The posterior mediation effect size.}
   \item{Post_Med_Prob}{The posterior Pr(Mediation) of this SNP.}
   \item{Post_Med_Prob2}{If `Post_Med_Prob` is less than `Pr_Med_cut`, set `Post_Med_Prob2` to 0. Otherwise, it's the same as `Post_Med_Prob`.}
   \item{SD_Post_CS_Prob}{Standardized Pr(CS) where Pr(CS) = Pr(Model)*Pr(Mediation)}
   \item{CumSum_Porb}{The cumulative `SD_Post_CS_Prob`. Note that the table is ordered by descending `SD_Post_CS_Prob`.}
   \item{EmpiricalCut}{The empirical coverage of this CS (should be >= requested `coverage`).}
   \item{CS_in}{A logical variable indicating whether this CS_SNP is included in this CS or not.}
   \item{index_SNP}{The name of the index SNP.}
}
}
\description{
Construct mJAM credible set based for selected index SNP
}
\author{
Jiayi Shen
}
