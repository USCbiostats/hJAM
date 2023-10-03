
#' Run mJAM with SuSiE
#'
#' @description fitting mJAM-SuSiE
#'
#' @param Marg_Result A data frame with marginal summary statistics from all studies. Col1: SNP name; Col2: Effect sizes from study #1; Col3: Std Errors of effect sizes from study #1; ...
#' @param EAF_Result A data frame with effect allele frequency (EAF) from all studies. Col1: SNP name; Col2: EAF from study #1; Col3: EAF from study #2; ...
#' @param N_GWAS A vector of sample sizes in all original GWAS studies.
#' @param X_ref A list of matrices with individual-level SNP dosage data in each study/population.
#' @param filter_rare A logical variable indicating whether to filter rare SNPs before the analysis. Default is `FALSE.` If `TRUE`, then please specify `rare_freq`.
#' @param rare_freq A vector of frequencies between 0 and 0.5 to specify the minor allele frequency cut-off if you want to filter rare SNPs before the analysis. Please also set `filter_rare` to be TRUE. For example, if there are 3 populations, then rare_freq = c(0.01, 0, 0.01) means SNPs with MAF < 0.01 in pop 1 and MAF < 0.01 in pop 3 will be removed from analysis.
#' @param SuSiE_num_comp SuSiE argument. The maximum number of causal SNPs that you want to select. Default is 10.
#' @param SuSiE_coverage SuSiE argument. The required coverage of credible sets. Default is 0.95.
#' @param SuSiE_min_abs_corr SuSiE argument. Minimum absolute correlation allowed in a credible set.
#' @param max_iter SuSiE argument. Maximum iterations to perform.
#' @param estimate_residual_variance SuSiE argument. If `TRUE`, then the susie algorithm is updating residual variance estimate during iterations. If `FALSE`, then use the residual variance is a fixed value, which is usually var(Y).
#'
#' @importFrom utils install.packages installed.packages
#' @import tibble
#' @import dplyr
#'
#' @export
#'
#' @author Jiayi Shen
#'
#' @returns
#'
#' \describe{
#'    \item{summary}{A table of the SuSiE posterior inclusion probabilities (PIPs), posterior mean, and posterior sd of all SNPs.}
#'    \item{fit}{SuSiE fit object.}
#' }
#'


mJAM_SuSiE <- function(Marg_Result = NULL,
                       EAF_Result = NULL,
                       N_GWAS,
                       X_ref,
                       filter_rare = FALSE,
                       rare_freq = NULL,
                       SuSiE_num_comp = 10,
                       SuSiE_coverage = 0.95,
                       SuSiE_min_abs_corr = 0.5,
                       max_iter = 500,
                       estimate_residual_variance = F){

  ###############################################################################
  ## Additional check points are in progress... to-do list:
  ## 1. whether all populations have the same dimension & names of SNPs
  ## 2. whether the length of marginal.logOR, marginal.se,EAFs, p_cases, N_GWAS, X_ref matches

  # if("susieR" %in% rownames(installed.packages()) == FALSE) {install.packages("susieR")}
  # library(susieR); library(dplyr); library(tibble)
  ###############################################################################

  ## Set parameters
  N_SNP <- numSNPs <- numSNPs_wo_rare <- nrow(Marg_Result)

  ## Check whether columns of X_ref, EAF_Results and Marg_Result are aligned.
  if(!all(sapply(X_ref, function(x) identical(colnames(x), colnames(X_ref[[1]]))))){
    stop("Columns of X_ref are not aligned with each other. \n
        Please make sure SNPs are sorted in the exact order.")
  }
  if(!identical(Marg_Result$SNP, EAF_Result$SNP)){
    stop("Rows in Marg_Result are not aligned to EAF_Result in the exact order.")
  }
  if(!identical(Marg_Result$SNP, colnames(X_ref[[1]]))){
    stop("Rows in Marg_Result are not aligned to X_ref in the exact order.")
  }


  ## if filter_rare, then remove rare SNPs
  if(typeof(filter_rare)!="logical"){
    stop("Please specify filter_rare to be either TRUE or FALSE.")
  }
  if(filter_rare == TRUE){
    ## check rare frequency cutoff
    if(is.null(rare_freq)){
      stop("Please specify a vector of minor allele frequency thresholds (between 0 and 0.5).")
    }else if(length(rare_freq)!=length(X_ref)){
      stop("Length of rare_freq does not match with the number of populations")
    }else{
      rare_freq_matrix = matrix(rep(rare_freq,numSNPs), byrow = TRUE, nrow = numSNPs, ncol = length(rare_freq))
      ## filter SNPs whose summary stat MAF < rare_freq in "ALL" population
      rare_freq_true1 = abs(EAF_Result[,2:ncol(EAF_Result)]-0.5)>0.5-rare_freq_matrix
      rare_filter_id1 = which(rowMeans(rare_freq_true1,na.rm=T)==1)
      ## filter SNPs whose reference MAF < rare_freq in "ALL" population
      reference_EAF = matrix(NA, nrow = numSNPs, ncol = length(X_ref))
      for(k in 1:length(X_ref)){
        reference_EAF[,k] <- colMeans(X_ref[[k]],na.rm=T)/2
      }
      rare_freq_true2 = abs(reference_EAF-0.5)>0.5-rare_freq_matrix
      rare_filter_id2 = which(rowMeans(rare_freq_true2,na.rm=T)==1)
      rare_filter_id = union(rare_filter_id1, rare_filter_id2)
      if(length(rare_filter_id)>0){
        ## filter Marg_Result and EAF_Result and X_ref
        Marg_Result = Marg_Result[-rare_filter_id,]
        EAF_Result = EAF_Result[-rare_filter_id,]
        for(i in 1:length(X_ref)){X_ref[[i]] = X_ref[[i]][,-rare_filter_id]}
        numSNPs_wo_rare = numSNPs - length(rare_filter_id)
      }
    }
  }

  Original_Input_Dosage <- X_ref
  Input_MarglogOR <- Input_MargSElogOR <- Input_MAF <- feMeta_w <- feMeta_w_x_beta <- vector("list", length(X_ref))
  for(pop in 1:length(X_ref)){
    Input_MarglogOR[[pop]] <- Marg_Result[,2*pop]
    Input_MargSElogOR[[pop]] <- Marg_Result[,2*pop+1]
    Input_MAF[[pop]] <- EAF_Result[,1+pop]
    names(Input_MarglogOR[[pop]]) <- Marg_Result[,1]
    names(Input_MargSElogOR[[pop]]) <- Marg_Result[,1]
    names(Input_MAF[[pop]]) <- EAF_Result[,1]
    feMeta_w[[pop]] <- 1/Marg_Result[,2*pop+1]^2
    feMeta_w_x_beta[[pop]] <- feMeta_w[[pop]]*Marg_Result[,2*pop]
  }
  feMeta_se <- sqrt(1/rowSums(do.call(cbind, feMeta_w), na.rm = TRUE))
  feMeta_mean <- rowSums(do.call(cbind, feMeta_w_x_beta), na.rm = TRUE)/rowSums(do.call(cbind, feMeta_w), na.rm = TRUE)


  ## --- Check missing data
  numEthnic <- length(Input_MarglogOR)
  missing_tbl <- tibble(missing_ethnic_idx = integer(), missing_snp_idx = integer())

  for(i in 1:numEthnic){
    temp_missing_gwas <- union(which(is.na(Input_MarglogOR[[i]])), which(is.na(Input_MargSElogOR[[i]])))
    temp_missing_gwas <- union(temp_missing_gwas, which(is.na(Input_MAF[[i]])))
    temp_missing_dosage <- which(colSums(is.na(X_ref[[i]]))>0)
    if(length(temp_missing_gwas)>0){
      missing_tbl <- missing_tbl %>% add_row(missing_ethnic_idx =i, missing_snp_idx = temp_missing_gwas)
    }
    if(length(temp_missing_dosage)>0){
      missing_tbl <- missing_tbl %>% add_row(missing_ethnic_idx =i, missing_snp_idx = temp_missing_dosage)
    }
  }
  missing_tbl <- distinct(missing_tbl)

  ## --- Transform logOR to betas
  # if (logORToBeta){
  #   linear.effects <- vector("list", numEthnic)
  #   for(i in 1:numEthnic){
  #     # Convert logOR to linear effects
  #     linear.effects[[i]] <- LogisticToLinearEffects(log.ors = marginal.betas[[i]],
  #                                                    log.or.ses = marginal.se[[i]],
  #                                                    snp.genotype.sds = NULL,
  #                                                    mafs = EAFs[[i]],
  #                                                    n = N_GWAS[i],
  #                                                    p.cases = p_cases[i])
  #     marginal.betas[[i]] <- linear.effects[[i]]$linear.beta.hats
  #     marginal.se[[i]] <- linear.effects[[i]]$linear.beta.ses
  #   }
  # }

  ## -- Calculate sufficient statistics
  GItGI <- GIty <- yty <- yty_med <- vector("list", numEthnic)
  for (i in 1:numEthnic){
    if (nrow(missing_tbl) > 0 && i%in%missing_tbl$missing_ethnic_idx){
      ## --- Get missing SNP index
      temp_missing_snp_idx <- filter(missing_tbl, missing_ethnic_idx == i) %>% pull(missing_snp_idx)
      ## --- Get GItGI, GIty, yty with complete SNP
      temp_X_ref <- X_ref[[i]][,-temp_missing_snp_idx]
      temp_MAFs <- Input_MAF[[i]][-temp_missing_snp_idx]
      temp_GItGI <- get_XtX(N_outcome = N_GWAS[i], Gl = temp_X_ref, maf = temp_MAFs)
      ##
      temp.marginal.betas <- Input_MarglogOR[[i]][-temp_missing_snp_idx]
      temp.marginal.se <- Input_MargSElogOR[[i]][-temp_missing_snp_idx]
      temp.GIty <- get_z(maf = temp_MAFs, betas = temp.marginal.betas, N_outcome = N_GWAS[i])
      ##
      # yty[[i]] <- get_yty(maf = temp_MAFs, N_outcome = N_GWAS[i], betas = temp.marginal.betas, betas.se = temp.marginal.se)
      ## --- Fill in missing SNPs with zeros
      GItGI[[i]] <- matrix(0, nrow = numSNPs_wo_rare, ncol = numSNPs_wo_rare)
      GItGI[[i]][-temp_missing_snp_idx, -temp_missing_snp_idx] <- temp_GItGI
      GIty[[i]] <- rep(0, numSNPs_wo_rare)
      GIty[[i]][-temp_missing_snp_idx] <- temp.GIty
    }else{
      GItGI[[i]] <- get_XtX(N_outcome = N_GWAS[i], Gl = X_ref[[i]],
                            maf = Input_MAF[[i]])
      GIty[[i]] <- get_z(maf = Input_MAF[[i]],
                         betas = Input_MarglogOR[[i]], N_outcome = N_GWAS[i])
    }
    temp_yty <- get_yty(maf = Input_MAF[[i]], N_outcome = N_GWAS[i],
                        betas = Input_MarglogOR[[i]],
                        betas.se = Input_MargSElogOR[[i]])
    yty[[i]] <- temp_yty$yTy.all
    yty_med[[i]] <- temp_yty$yTy.est
    colnames(GItGI[[i]]) <- rownames(GItGI[[i]]) <- Marg_Result[,1]
    names(GIty[[i]]) <- names(yty[[i]]) <- Marg_Result[,1]
  }

  susie_in_XtX <- Reduce("+", GItGI)
  susie_in_Xty <- Reduce("+", GIty)
  susie_in_yty <- Reduce("+", yty_med)

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

  pip <- susieR::susie_get_pip(susie_fit)
  post_mean <- susie_get_posterior_mean_v2(susie_fit)
  post_sd <- susie_get_posterior_sd_v2(susie_fit)
  summary_table <- data.frame(pip = pip,post_mean = post_mean,post_sd = post_sd,stringsAsFactors = FALSE)
  rownames(summary_table) <- Marg_Result[,1]

  return(list(summary = summary_table, fit = susie_fit))
}
