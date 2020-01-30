#' Example reference data of hJAM
#'
#' The real data example from hJAM paper
#'
#' @docType data
#'
#' @format The \code{Gl} object is a data matrix with 2467 individual of 210 SNPs from 1000 Genome project.
#' @references Consortium GP. A global reference for human genetic variation. Nature 2015; 526: 68.
"Gl"

#' Example beta list of hJAM
#'
#' @docType data
#' @format The \code{betas.Gy} is the beta vector in the hJAM model: the association estimates between 210 SNPs and myocardial infarction. The summary data was collected from UK Biobank (n=459,324).
#' @references Sudlow C, Gallacher J, Allen N, et al. UK biobank: an open access resource for identifying the causes of a wide range of complex diseases of middle and old age. PLoS Med 2015; 12: e1001779.
"betas.Gy"

#' Example conditional A matrix of hJAM
#'
#' @format The \code{conditional_A} is the conditional estimates alpha matrix in the hJAM model: the association estimates between 210 SNPs and body mass index (BMI) and type 2 diabetes (T2D). The summary data was collected from GIANT consortium (n=339,224) and DIAGRAM+GERA+UKB (n=659316) for BMI and T2D, respectively. We converted it from marginal_A, using \code{get_cond_A} function in hJAM package.
#' @references 1. Locke AE, Kahali B, Berndt SI, et al. Genetic studies of body mass index yield new insights for obesity biology. Nature 2015; 518: 197-206. 2. Xue A, Wu Y, Zhu Z, et al. Genome-wide association analyses identify 143 risk variants and putative regulatory mechanisms for type 2 diabetes. Nat Commun 2018; 9: 2941.
"conditional_A"

#' Example marginal A matrix of hJAM
#'
#' @format The \code{marginal_A} is the marginal estimates alpha matrix in the hJAM model: the association estimates between 210 SNPs and body mass index (BMI) and type 2 diabetes (T2D). The summary data was collected from GIANT consortium (n=339,224) and DIAGRAM+GERA+UKB (n=659316) for BMI and T2D, respectively.
#' @references 1. Locke AE, Kahali B, Berndt SI, et al. Genetic studies of body mass index yield new insights for obesity biology. Nature 2015; 518: 197-206. 2. Xue A, Wu Y, Zhu Z, et al. Genome-wide association analyses identify 143 risk variants and putative regulatory mechanisms for type 2 diabetes. Nat Commun 2018; 9: 2941.
"marginal_A"

#' Example SNPs' information of hJAM
#'
#' @format The \code{SNPs_info} is the information of the 210 SNPs that we used in this data example. It includes three columns: the rsID, major allele, and minor allele frequency of each SNP. The minor allele frequencies were calculated in the 503 European-ancestry subjects in 1000 Genome project.
#' @references Consortium GP. A global reference for human genetic variation. Nature 2015; 526: 68.
"SNPs_info"
