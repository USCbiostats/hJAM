### Manual LD pruning function 

## Input: 
## 1. a list of correlation matrix of all SNPs (GItGI_LA_R)
## 2. a list of MAF of all SNPs 
## 3. Target SNP ID and the IDs of Testing SNPs 

# target = target_id
# testing = testing_ids
# R <- list(GItGI_EUR_R, GItGI_AA_R, GItGI_LA_R, GItGI_Asian_R)
# within_thre = 0.95
# across_thre = 0.80

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

