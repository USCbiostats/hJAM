#' Example data of hJAM
#'
#' Real data for BMI/T2D on the risk of Myocardial infarction
#'
#' @format The \code{MI} object is a set of data sets which was used to estimate the causal effect of body mass index and type 2 diabetes on the risk of myocardial infarction.
#' \describe{
#'  \item{MI.marginal.Amatrix}{The marginal \eqn{\hat{A}} matrix. Column one and two are the marginal estimates of the SNPs on body mass index from GIANT consortium (n = 339,224) (Locke et al., 2015) and type 2 diabetes from DIAGRAM+GERA+UKB (n = 659,316) (Xue et al., 2018), respectively}
#'  \item{MI.Amatrix}{The conditional \eqn{\hat{A}} matrix composed by JAM and the marginal \eqn{\hat{A}} matrix. Column one and two are the conditional effect estimates of the SNPs on body mass index and type 2 diabetes, respectively.}
#'  \item{MI.Geno}{The reference genotype data from the European-ancestry population in 1000 Genome Project (Consortium, 2015).}
#'  \item{MI.betas.gwas}{The b vector. The association estimates between selected SNPs and the risk of myocardial infarction from UK Biobank (Sudlow et al., 2015).}
#'  \item{MI.SNPs_info}{The SNP information. Five columns included: the RSID, reference allele, reference allele frequency, if BMI significant and if T2D significant. The last two columns are indicator variables for the SNPs which are genome-wide significant associated with BMI/T2D.}
#' }
#' @references Consortium GP. A global reference for human genetic variation. Nature 2015; 526: 68.
#' @references Locke, Adam E., et al. Genetic studies of body mass index yield new insights for obesity biology. Nature 518.7538 (2015): 197-206.
#' @references Xue, Angli, et al. Genome-wide association analyses identify 143 risk variants and putative regulatory mechanisms for type 2 diabetes. Nature communications 9.1 (2018): 1-14.
#' @references Sudlow, Cathie, et al. UK biobank: an open access resource for identifying the causes of a wide range of complex diseases of middle and old age. Plos med 12.3 (2015): e1001779.
#' @name MI
#' @aliases MI.marginal.Amatrix MI.Amatrix MI.Geno MI.betas.gwas MI.SNPs_info
NULL

#' Real data for selecting the genes on chromosome 10 for the risk of prostate cancer
#'
#' @format The \code{GTEx.PrCa} is a set of data sets which was applied for selecting the genes on chromosome 10 for the risk of prostate cancer
#' \describe{
#'  \item{GTEx.PrCa.marginal.Amatrix}{The marginal \eqn{\hat{A}} matrix with 158 genes and 182 eQTLs. The raw data was downloaded from GTEx analysis v7 (https://gtexportal.org/home/datasets). Priority Pruner was used to select the independent eQTLs. We used this matrix for MR-BMA implementation.)}
#'  \item{GTEx.PrCa.Amatrix}{The conditional \eqn{\hat{A}} matrix with 167 genes and 447 eQTLs, which was composed by SuSiE JAM and the raw data of \eqn{\hat{A}} matrix.}
#'  \item{GTEx.PrCa.Geno}{The reference genotype data for the 447 eQTLs from the European-ancestry population in 1000 Genome Project (Consortium, 2015)}
#'  \item{GTEx.PrCa.betas.gwas}{The b vector. The association estimates between eQTLs and the risk of prostate cancer from (Schumacher et al., 2018)}
#'  \item{GTEx.PrCa.betas.se.gwas}{The se(b) vector from (Schumacher et al., 2018)}
#'  \item{GTEx.PrCa.pvalue.gwas}{The pvalues vector of the association estimates between selected SNPs and the risk of prostate cancer from (Schumacher et al., 2018)}
#'  \item{GTEx.PrCa.maf.gwas}{The vector of the effect allele frequency of the SNPs from (Schumacher et al., 2018)}
#' }
#' @references Consortium GP. A global reference for human genetic variation. Nature 2015; 526: 68.
#' @references Lonsdale, John, et al. The genotype-tissue expression (GTEx) project. Nature genetics 45.6 (2013): 580-585.
#' @references Schumacher, Fredrick R., et al. Association analyses of more than 140,000 men identify 63 new prostate cancer susceptibility loci. Nature genetics 50.7 (2018): 928-936.
#' @name GTEx.PrCa
#' @aliases GTEx.PrCa.marginal.Amatrix GTEx.PrCa.Amatrix GTEx.PrCa.Geno GTEx.PrCa.betas.gwas GTEx.PrCa.betas.se.gwas GTEx.PrCa.pvalues.gwas GTEx.PrCa.maf.gwas GTEx.PrCa.rsid GTEx.PrCa.marginal.selected
NULL

#' Real data for selecting the metabolites for the risk of prostate cancer

#' @format The \code{PrCa.lipids} is a set of data sets which was for selecting the metabolites for the risk of prostate cancer
#' \describe{
#'  \item{PrCa.lipids.marginal.Amatrix}{The marginal \eqn{\hat{A}} matrix with 118 metabolites and 144 SNPs. This data is directly adapted from https://github.com/verena-zuber/demo_AMD (Zuber et al., 2020)}
#'  \item{PrCa.lipids.Amatrix}{The conditional \eqn{\hat{A}} matrix with 118 metabolites and 144 SNPs, which was composed by SuSiE JAM and the marginal \eqn{\hat{A}} matrix.}
#'  \item{PrCa.lipids.Geno}{The reference genotype data for the 144 SNPs from the European-ancestry population in 1000 Genome Project (Consortium, 2015).}
#'  \item{PrCa.lipids.betas.gwas}{The b vector. The association estimates between selected SNPs and the risk of prostate cancer from (Schumacher et al., 2018)}
#'  \item{PrCa.lipids.betas.se.gwas}{The se(b) vector from (Schumacher et al., 2018)}
#'  \item{PrCa.lipids.pvalue.gwas}{The pvalues vector of the association estimates between selected SNPs and the risk of prostate cancer from (Schumacher et al., 2018)}
#'  \item{PrCa.lipids.maf.gwas}{The vector of the effect allele frequency of the SNPs from (Schumacher et al., 2018)}
#'  \item{PrCa.lipids.rsid}{The RSID of the SNPs.}
#' }
#' @references Consortium GP. A global reference for human genetic variation. Nature 2015; 526: 68.
#' @references Zuber, Verena, et al. Selecting likely causal risk factors from high-throughput experiments using multivariable Mendelian randomization. Nature communications 11.1 (2020): 1-11.
#' @references Schumacher, Fredrick R., et al. Association analyses of more than 140,000 men identify 63 new prostate cancer susceptibility loci. Nature genetics 50.7 (2018): 928-936.
#' @name PrCa.lipids
#' @aliases PrCa.lipids.marginal.Amatrix PrCa.lipids.Amatrix PrCa.lipids.Geno PrCa.lipids.betas.gwas PrCa.lipids.betas.se.gwas PrCa.lipids.pvalue.gwas PrCa.lipids.maf.gwas PrCa.lipids.rsid
NULL
