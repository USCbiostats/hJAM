
#' Pruning SNPs based on LD
#'
#'
#' @param target Target SNP ID.
#' @param testing IDs of SNPs to be tested.
#' @param R a list of correlation matrix of all SNPs.
#' @param within_thre threshold of r2 with selected index SNP(s) within a single population. If a SNP's correlation with any selected index SNP is greater than this threshold in at least one population, it will be excluded from subsequent rounds of index SNP selection.
#' @param across_thre threshold of r2 with selected index SNP(s) across all populations. If a SNP's correlation with any selected index SNP is greater than this threshold in all populations, it will be excluded from subsequent rounds of index SNP selection.
#'
#' @author Jiayi Shen
#'
#' @returns
#'
#' \describe{
#'    \item{remove_within}{SNP IDs to be pruned due to high within-population correlation}
#'    \item{remove_across}{SNP IDs to be pruned due to high across-population correlation}
#' }
#'


mJAM_LDpruning = function(target, testing, R, within_thre = 0.95, across_thre = 0.80){

  numEthnic <- length(R)
  high_corr_within <- high_corr_across <- vector(mode = "list", length = numEthnic)

  for(i in 1:numEthnic){
    high_corr_within[[i]] <- which(abs(R[[i]][target,]) > within_thre)
    high_corr_within[[i]] <- intersect(testing, high_corr_within[[i]])
    high_corr_across[[i]] <- which(abs(R[[i]][target,]) > across_thre)
    high_corr_across[[i]] <- intersect(testing, high_corr_across[[i]])
  }

  remove_within <- Reduce(union, high_corr_within)
  remove_across <- Reduce(intersect, high_corr_across)

  return(list(remove_within = remove_within, remove_across = remove_across))
}

