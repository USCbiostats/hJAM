###############################################################################
## This script is for mJAM (multi-ethnic Joint Analysis of Marginal statistics) 
## with SuSiE (sufficient statistic version) as the model selection method.
## 
## Section 1: Functions to assist intermediate steps in the main mJAM-SuSiE function
## Section 2: The main function of mJAM-SuSiE
## Section 3: Functions to draw inferences from a fitted mJAM-SuSiE model. 
## 

##########################  Section 1  #################################
## Function 1.1 -- Get transformed statistics: z, or Xty
get_z <- function(maf, betas, N_outcome){
  
  ## follows JAM supplementary material section 1
  n0 = N_outcome*(1-maf)^2
  n1 = N_outcome*2*maf*(1-maf)
  n2 = N_outcome*maf^2
  
  y0 = -(n1*betas+2*n2*betas)/(n0+n1+n2)
  y1 = y0+betas
  y2 = y0+2*betas
  z = n1*y1 + 2*n2*y2
  
  return(z)
}


get_Xty <- function(maf, betas, N_outcome){
  
  Dj <- 2*maf*(1-maf)*N_outcome
  Xty <- diag(Dj) %*% betas
  
  return(Xty)
}

## Function 1.2 -- Get var-cov matrix: XtX
get_XtX <- function(N_outcome, Gl, maf){
  
  ## --- center Gl to have mean 0
  G0 <- scale(Gl, center = TRUE, scale = FALSE)
  G0_t_G0 <- t(G0)%*%G0
  
  ## --- Modify G0'G0 if the sample sizes of Gl and Gx are different
  if (nrow(G0_t_G0)>1){
    Dj <- 2*maf*(1-maf)*N_outcome
    D_sqrt <- diag(sqrt(Dj))
    Dw_sqrt_inv <- diag(1/sqrt(diag(G0_t_G0)))
    G0_t_G0.scaled <- D_sqrt %*% Dw_sqrt_inv  %*% G0_t_G0 %*% Dw_sqrt_inv %*% D_sqrt
    ## Add a ridge term in case G0'G0 is singular
    if(matrixcalc::is.singular.matrix(G0_t_G0.scaled)){
      ridgeValue <-  min(1,min(diag(G0_t_G0.scaled)*.001))
      G0_t_G0.scaled <-  G0_t_G0.scaled + ridgeValue*diag(dim(G0_t_G0.scaled)[2])
      message("G0_t_G0.scaled is singular. Ridge term added.")
    }
  }else{
    Dj <- 2*maf*(1-maf)*N_outcome
    D_sqrt <- sqrt(Dj)
    Dw_sqrt_inv <- 1/sqrt(G0_t_G0)
    G0_t_G0.scaled <- D_sqrt * Dw_sqrt_inv * G0_t_G0 * Dw_sqrt_inv * D_sqrt
  }
  
  return(G0_t_G0.scaled)
}


## Function 1.3 --- Obtain sum of squares, yty
get_yty <- function(maf, N_outcome, betas, betas.se){
  
  ## Follows Yang's Nature paper (2012)
  Dj <- 2*maf*(1-maf)*N_outcome
  Sj2 <- betas.se^2
  yTy.all <- Dj*Sj2*(N_outcome-1)+Dj*betas^2
  yTy.est <- median(Dj*Sj2*(N_outcome-1)+Dj*betas^2, na.rm = TRUE)
  # yTy.all[is.na(yTy.all)] <- median(Dj*Sj2*(N_outcome-1)+Dj*betas^2, na.rm = TRUE)
  yTy.all[is.na(yTy.all)] <-  0
  
  return(list(yTy.est = yTy.est, yTy.all = yTy.all))
}



## Function 1.4 --- Transform logistic effect to linear effect 
## (pasted from R2BGLiMS::JAM_LogisticToLinearEffects)
## Reference: Benner 2015, FINEMAP
LogisticToLinearEffects <- function(
  log.ors = NULL,
  log.or.ses = NULL,
  snp.genotype.sds = NULL,
  mafs = NULL,
  n = NULL,
  p.cases = NULL
) {
  
  # Standardised least squares estimate is signed z-score/sqrt(n)
  standardised.least.squares.effect <- (log.ors/log.or.ses)/sqrt(n)
  
  if (!is.null(mafs) & !is.null(snp.genotype.sds)) {
    cat("snp.sds and mafs were provided. snp.genotype.sds will be used as the preferred option.\n")
    mafs <- NULL
  }
  
  if (!is.null(mafs)) {
    # Multiply by trait SD for effect on trait scale and divide by SNP SD for per allele effect
    linear.beta.hats <- standardised.least.squares.effect*sqrt(p.cases*(1-p.cases))/sqrt(2*mafs*(1-mafs))
    linear.beta.ses <- sqrt(p.cases*(1-p.cases))/(sqrt(n)*sqrt(2*mafs*(1-mafs)))
  }
  
  if (!is.null(snp.genotype.sds)) {
    # Multiply by trait SD for effect on trait scale and divide by SNP SD for per allele effect
    linear.beta.hats <- standardised.least.squares.effect*sqrt(p.cases*(1-p.cases))/snp.genotype.sds
    linear.beta.ses <- sqrt(p.cases*(1-p.cases))/(sqrt(n)*snp.genotype.sds)
  }
  
  return(list(linear.beta.hats = linear.beta.hats, linear.beta.ses = linear.beta.ses))
}

##########################  Section 2  #################################
## -- This is the main function to run mJAM-SuSiE. 
mJAM_SuSiE <- function(marginal.betas = NULL, 
                       marginal.se = NULL, 
                       EAFs, 
                       N_GWAS,
                       X_ref,
                       SNP_names = NULL,
                       logORToBeta = FALSE,
                       p_cases = NULL, 
                       replace_missing = TRUE,
                       SuSiE_num_comp = 10, 
                       SuSiE_coverage = 0.95, 
                       SuSiE_min_abs_corr = 0.5, 
                       max_iter = 500,
                       estimate_residual_variance = F){
  
  ###############################################################################
  ## Additional check points are in progress... to-do list:
  ## 1. whether all populations have the same dimension & names of SNPs 
  ## 2. whether the length of marginal.logOR, marginal.se,EAFs, p_cases, N_GWAS, X_ref matches
  
  if("susieR" %in% rownames(installed.packages()) == FALSE) {install.packages("susieR")}
  if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
  library(susieR); library(tidyverse)
  ###############################################################################
  
  ## --- uncomment the following only for the convenience of debugging 
  # marginal.betas = list(MargBeta$MargBeta_P1, MargBeta$MargBeta_P2, MargBeta$MargBeta_P3)
  # marginal.se = list(MargBeta$MargSEBeta_P1, MargBeta$MargSEBeta_P2, MargBeta$MargSEBeta_P3)
  # EAFs = list(EAF$EAF_P1, EAF$EAF_P2, EAF$EAF_P3)
  # N_GWAS = c(5000, 5000, 5000)
  # X_ref = list(RefDosage_P1,RefDosage_P2,RefDosage_P3)
  # logORToBeta = FALSE
  # p_cases = c(0.5,0.5,0.5)
  # replace_missing = TRUE
  # SNP_names = MargBeta$SNP
  # SuSiE_num_comp = 10
  # SuSiE_coverage = 0.95
  # SuSiE_min_abs_corr = 0.5
  # max_iter = 500
  # estimate_residual_variance = F
  
  ## --- Check missing data 
  numEthnic <- length(marginal.betas)
  numSNPs <- length(marginal.betas[[1]])
  missing_tbl <- tibble(missing_ethnic_idx = integer(), missing_snp_idx = integer())
  
  for(i in 1:numEthnic){
    temp_missing_gwas <- union(which(is.na(marginal.betas[[i]])), which(is.na(marginal.se[[i]])))
    temp_missing_gwas <- union(temp_missing_gwas,which(is.na(EAFs[[i]])))
    temp_missing_dosage <- which(colSums(is.na(X_ref[[i]]))>0)
    if(length(temp_missing_gwas)>0){
      missing_tbl <- missing_tbl %>% add_row(missing_ethnic_idx =i, missing_snp_idx = temp_missing_gwas)
    }
    if(length(temp_missing_dosage)>0){
      missing_tbl <- missing_tbl %>% add_row(missing_ethnic_idx =i, missing_snp_idx = temp_missing_dosage)
    }
  }
  
  if(nrow(missing_tbl) > 0 && replace_missing == FALSE){
    stop("Missing values in at least one population. Please specify replace_missing.")
  }
  
  ## --- Transform logOR to betas
  if (logORToBeta){
    linear.effects <- vector("list", numEthnic)
    for(i in 1:numEthnic){
      # Convert logOR to linear effects 
      linear.effects[[i]] <- LogisticToLinearEffects(log.ors = marginal.betas[[i]],
                                                     log.or.ses = marginal.se[[i]],
                                                     snp.genotype.sds = NULL,
                                                     mafs = EAFs[[i]],
                                                     n = N_GWAS[i],
                                                     p.cases = p_cases[i])
      marginal.betas[[i]] <- linear.effects[[i]]$linear.beta.hats
      marginal.se[[i]] <- linear.effects[[i]]$linear.beta.ses
    }
  }
  
  ## --- Get the SuSiE inputs: I'G'GI,I'G'y and y'y
  GItGI <- GIty <- yty <- vector("list", numEthnic)
  
  if(nrow(missing_tbl) > 0 && replace_missing){
    ################################################################################
    for (i in 1:numEthnic){
      if (i%in%missing_tbl$missing_ethnic_idx){
        ## --- Get missing SNP index
        temp_missing_snp_idx <- filter(missing_tbl, missing_ethnic_idx == i) %>% pull(missing_snp_idx)
        ## --- Get GItGI, GIty, yty with complete SNP
        temp_X_ref <- X_ref[[i]][,-temp_missing_snp_idx]
        temp_EAFs <- EAFs[[i]][-temp_missing_snp_idx]
        temp_GItGI <- get_XtX(N_outcome = N_GWAS[i], Gl = temp_X_ref, maf = temp_EAFs)
        ## 
        temp.marginal.betas <- marginal.betas[[i]][-temp_missing_snp_idx]
        temp.marginal.se <- marginal.se[[i]][-temp_missing_snp_idx]
        temp.GIty <- get_z(maf = temp_EAFs, betas = temp.marginal.betas, N_outcome = N_GWAS[i])
        ##
        yty[[i]] <- get_yty(maf = temp_EAFs, N_outcome = N_GWAS[i], betas = temp.marginal.betas, betas.se = temp.marginal.se)$yTy.est
        ## --- Fill in missing SNPs with zeros
        GItGI[[i]] <- matrix(0, nrow = numSNPs, ncol = numSNPs)
        GItGI[[i]][-temp_missing_snp_idx, -temp_missing_snp_idx] <- temp_GItGI
        GIty[[i]] <- rep(0, numSNPs)
        GIty[[i]][-temp_missing_snp_idx] <- temp.GIty
      }else{
        GItGI[[i]] <- get_XtX(N_outcome = N_GWAS[i], Gl = X_ref[[i]], maf = EAFs[[i]])
        GIty[[i]] <- get_z(maf = EAFs[[i]], betas = marginal.betas[[i]], N_outcome = N_GWAS[i])
        yty[[i]] <- get_yty(maf = EAFs[[i]], N_outcome = N_GWAS[i], betas = marginal.betas[[i]], betas.se = marginal.se[[i]])$yTy.est
      }
    }
    ################################################################################
  }else{
    for (i in 1:numEthnic){
      GItGI[[i]] <- get_XtX(N_outcome = N_GWAS[i], Gl = X_ref[[i]], maf = EAFs[[i]])
      GIty[[i]] <- get_z(maf = EAFs[[i]], betas = marginal.betas[[i]], N_outcome = N_GWAS[i])
      yty[[i]] <- get_yty(maf = EAFs[[i]], N_outcome = N_GWAS[i], betas = marginal.betas[[i]], betas.se = marginal.se[[i]])$yTy.est
    }
  }
  

  susie_in_XtX <- Reduce("+", GItGI)
  susie_in_Xty <- Reduce("+", GIty)
  susie_in_yty <- Reduce("+", yty)
  
  ## --- Run SuSiE sufficient statistics version
  susie_fit <- susieR::susie_suff_stat(XtX = susie_in_XtX, 
                                       Xty = susie_in_Xty,
                                       n = sum(N_GWAS), 
                                       yty = susie_in_yty, 
                                       L = SuSiE_num_comp,
                                       min_abs_corr = SuSiE_min_abs_corr, 
                                       max_iter = max_iter, 
                                       coverage= SuSiE_coverage,
                                       standardize = T, 
                                       estimate_residual_variance = estimate_residual_variance, 
                                       verbose = F)
  
  ## --- Outputs
  # 1. posterior probability of each SNP
  # 2. 95% credible sets 
  # 3. (optional) posterior mean effect of each SNP 
  
  pip <- susie_get_pip(susie_fit)
  post_mean <- susie_get_posterior_mean_v2(susie_fit)
  post_sd <- susie_get_posterior_sd_v2(susie_fit)
  summary_table <- data.frame(pip = pip,post_mean = post_mean,post_sd = post_sd,stringsAsFactors = FALSE)
  rownames(summary_table) <- SNP_names
  
  return(list(summary = summary_table, fit = susie_fit))
}


##########################  Section 3  #################################
## Function 3.1 -- Calculate posterior mean from SuSiE output 
##   (Modified from susieR::susie_get_posterior_mean; 
##   v2 fixed the error related to colSums(). )
susie_get_posterior_mean_v2 <- function (res, prior_tol = 1e-09) 
{
  if (is.numeric(res$V)) {
    include_idx = which(res$V > prior_tol)
  }
  else {include_idx = 1:nrow(res$alpha)}
  
  if (length(include_idx) > 1) {
    return(colSums((res$alpha * res$mu)[include_idx, ])/res$X_column_scale_factors)
  } else if (length(include_idx) == 1) {
    return((res$alpha * res$mu)[include_idx, ]/res$X_column_scale_factors)
  }else {
    return(numeric(ncol(res$mu)))
  } 
}

## Function 3.2 --- Calculate posterior standard deviation from SuSiE output 
##     (Modified from susieR::susie_get_posterior_sd; 
##      v2 fixed the error related to colSums(). )
susie_get_posterior_sd_v2 <- function (res, prior_tol = 1e-09) 
{
  if (is.numeric(res$V)) 
    include_idx = which(res$V > prior_tol)
  else include_idx = 1:nrow(res$alpha)
  if (length(include_idx) > 1) {
    return(sqrt(colSums((res$alpha * res$mu2 - (res$alpha * 
                                                  res$mu)^2)[include_idx, ]))/(res$X_column_scale_factors))
  } else if (length(include_idx) == 1){
    return(sqrt((res$alpha * res$mu2 - (res$alpha * res$mu)^2)[include_idx, ])/(res$X_column_scale_factors))
  } else {return(numeric(ncol(res$mu)))}
}

## Function 3.3 --- Tidy and present the credible set from a mJAM-SuSiE output. 
mJAM_SuSiE_get_cs <- function(susie_fit, SNP_names, 
                              coverage = 0.95){
  
  # cs_output <- susieR::susie_get_cs(susie_fit, coverage = coverage)
  cs_output <- susie_fit$sets
  if (is.null(cs_output$cs)){
    cs_summary <- data.frame(index = NA, 
                             coverage = coverage,
                             CS_size = 0,
                             index_SNP = "None Selected",
                             CS_SNP = "None Selected"
                             )
  }else{
    cs_summary <- data.frame(index = character(),coverage = double(),
                             CS_size = integer(),index_SNP_id = integer(), 
                             CS_SNP_id = integer())
    for(r in 1:length(cs_output$cs)){
      temp_cs_summary <- data.frame(index = names(cs_output$cs)[r],
                                    coverage = cs_output$coverage[r],
                                    CS_size = length(cs_output$cs[[r]]),
                                    index_SNP_id = cs_output$cs_index[r],
                                    CS_SNP_id = cs_output$cs[[r]]
                                    )
      cs_summary <- rbind(cs_summary,temp_cs_summary)
    }
    cs_summary <- cs_summary %>% 
      left_join(data.frame(index_SNP_id = 1:length(SNP_names), index_SNP = SNP_names), by = "index_SNP_id") %>% 
      left_join(data.frame(CS_SNP_id = 1:length(SNP_names), CS_SNP = SNP_names), by = "CS_SNP_id") %>% 
      dplyr::select(-c(index_SNP_id,CS_SNP_id))
      
  }

  return(cs_summary)
  
}



##########################  Section 4  #################################
## Function 4.1 -- Calculate null-based Bayes Factor using g-prior 
 
mJAM_get_BF <- function(marginal.betas, marginal.se,MAFs, X_ref,N_GWAS, selected_id, g = NULL){
  
  ## --- Check missing data
  numEthnic <- length(marginal.betas)
  numSNPs <- length(marginal.betas[[1]])
  missing_snp_idx <- missing_ethnic_idx <- NULL
  missing_tbl <- tibble(missing_ethnic_idx, missing_snp_idx)
  
  for(i in 1:numEthnic){
    temp_missing_gwas <- union(which(is.na(marginal.betas[[i]])), which(is.na(marginal.se[[i]])))
    temp_missing_dosage <- which(colSums(is.na(X_ref[[i]]))>0)
    if(length(temp_missing_gwas)>0){
      missing_tbl <- missing_tbl %>% add_row(missing_ethnic_idx =i, missing_snp_idx = temp_missing_gwas)
    }
    if(length(temp_missing_dosage)>0){
      missing_tbl <- missing_tbl %>% add_row(missing_ethnic_idx =i, missing_snp_idx = temp_missing_dosage)
    }
  }
  
  ## -- Calculate model-specific sufficient statistics 
  GItGI <- GIty <- yty <- vector("list", numEthnic)
  
  ## --- Specify which SNPs to include
  R2 <- rep(NA, numSNPs)
  All_SNPs <- 1:numSNPs
  in_id_list <- All_SNPs[!All_SNPs %in% selected_id]
  
  for (in_id in in_id_list){
    for (i in 1:numEthnic){
      if (nrow(missing_tbl) > 0 && i%in%missing_tbl$missing_ethnic_idx){
        ### TO-DO: Modify this part to handle missing SNP situations 
        # ## --- Get missing SNP index
        # temp_missing_snp_idx <- filter(missing_tbl, missing_ethnic_idx == i) %>% pull(missing_snp_idx)
        # ## --- Get GItGI, GIty, yty with complete SNP
        # temp_X_ref <- X_ref[[i]][,-temp_missing_snp_idx]
        # temp_MAFs <- MAFs[[i]][-temp_missing_snp_idx]
        # temp_GItGI <- get_XtX(N_outcome = N_GWAS[i], Gl = temp_X_ref, maf = temp_MAFs)
        # ##
        # temp.marginal.betas <- marginal.betas[[i]][-temp_missing_snp_idx]
        # temp.marginal.se <- marginal.se[[i]][-temp_missing_snp_idx]
        # temp.GIty <- get_z(maf = temp_MAFs, betas = temp.marginal.betas, N_outcome = N_GWAS[i])
        # ##
        # yty[[i]] <- get_yty(maf = temp_MAFs, N_outcome = N_GWAS[i], betas = temp.marginal.betas, betas.se = temp.marginal.se)
        # ## --- Fill in missing SNPs with zeros
        # GItGI[[i]] <- matrix(0, nrow = numSNPs, ncol = numSNPs)
        # GItGI[[i]][-temp_missing_snp_idx, -temp_missing_snp_idx] <- temp_GItGI
        # GIty[[i]] <- rep(0, numSNPs)
        # GIty[[i]][-temp_missing_snp_idx] <- temp.GIty
      }else{
        GItGI[[i]] <- get_XtX(N_outcome = N_GWAS[i], Gl = X_ref[[i]][,c(selected_id, in_id)], maf = MAFs[[i]][c(selected_id, in_id)])
        GIty[[i]] <- get_z(maf = MAFs[[i]][c(selected_id, in_id)], betas = marginal.betas[[i]][c(selected_id, in_id)], N_outcome = N_GWAS[i])
        yty[[i]] <- get_yty(maf = MAFs[[i]][c(selected_id, in_id)], N_outcome = N_GWAS[i], 
                            betas = marginal.betas[[i]][c(selected_id, in_id)], betas.se = marginal.se[[i]][c(selected_id, in_id)])
      }
    }
    
    susie_in_XtX <- Reduce("+", GItGI)
    susie_in_Xty <- Reduce("+", GIty)
    susie_in_yty <- Reduce("+", yty)
    
    ## --- Calculate multi-ethnic R-squared
    if(nrow(susie_in_XtX) == 1){
      b_joint <- susie_in_Xty / susie_in_XtX
      R2_est <- b_joint*susie_in_XtX*b_joint/susie_in_yty
      R2[in_id] <- ifelse(R2_est > 1, 1, R2_est)
    }else{
      D <- diag(diag(susie_in_XtX))
      b_marg <- solve(D) %*% susie_in_Xty
      b_joint <- solve(susie_in_XtX) %*% D %*% b_marg
      R2_est <- t(b_joint) %*% D %*% b_marg / susie_in_yty
      R2[in_id] <- ifelse(R2_est > 1, 1, R2_est)
    }
  }
  

  ## Specify the prior information 
  if(is.null(g)){g <- sum(N_GWAS)}
  
  ## Calculate BF
  BF <- vector(mode = "list", length = numSNPs)
  BF_max <- -Inf; which_BF_max <- NA
  p_model <- length(selected_id)+1
  
  for (in_id in in_id_list){
    log1 <- Rmpfr::mpfr(log(1+g), 50)
    log2 <- Rmpfr::mpfr(log(1+g*(1-R2[in_id])), 50)
    BF[[in_id]] <- 0.5*(sum(N_GWAS)-p_model-1)*log1 - 0.5*(sum(N_GWAS)-1)*log2
    if(BF[[in_id]]>BF_max){
      BF_max <- BF[[in_id]]
      which_BF_max <- in_id
    }
  }
  
  return(list(which_BF_max = which_BF_max, BF_max = BF_max, BF = BF, r2_max = R2[which_BF_max]))
}



## Function 4.3: Calculate within-CS/between-CS correlation 
## Input: Dosage matrix or correlation matrix (multiple population), SNPs within the CS
##        If CS1_SNPs = CS2_SNPs, then it calculates within-CS correlation. 
## Output: A data frame with empirical correlation between any pair of SNPs within the CS

mJAM_get_bwCS_corr <- function(X_ref, CS1_SNPs, CS2_SNPs){
  
  # CS1_SNPs = c(9,2,1); CS2_SNPs = c(34,2,32,36,1,40)
  N_Study <- length(X_ref)
  corr_VarName <- paste0("temp_corr_p", 1:N_Study)
  
  corr_df <- cbind(expand.grid(CS1_SNPs, CS2_SNPs), 
                   data.frame(matrix(NA, ncol = N_Study, nrow = length(CS2_SNPs)*length(CS2_SNPs))))
  
  for(j in 1:nrow(corr_df)){
    var1 <- corr_df$Var1[j]
    var2 <- corr_df$Var2[j]
    for (e in 1:N_Study){
      corr_df[j,2+e] <- cor(X_ref[[e]][,var1], X_ref[[e]][,var2])
    }
  }
  colnames(corr_df) <- c("SNP1","SNP2",paste0("P",1:N_Study, "_corr"))
  
  return(corr_df)
}
