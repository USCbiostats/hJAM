#' Run mJAM with Forward Selection
#'
#' @description fitting mJAM-Forward
#'
#' @param N_GWAS A vector of sample sizes in all original GWAS studies.
#' @param X_ref A list of matrices with individual-level SNP dosage data in each study/population. Each column corresponds to a SNP. Note that the columns name should match exactly to the SNP column in `Marg_Result` and `EAF_Result`. If certain SNP(s) is missing in dosage, then insert NAs in corresponding column(s).
#' @param Marg_Result A data frame with marginal summary statistics from all studies. Col1: SNP name; Col2: Effect sizes from study #1; Col3: Std Errors of effect sizes from study #1; ...
#' @param EAF_Result A data frame with effect allele frequency (EAF) from all studies. Col1: SNP name; Col2: EAF from study #1; Col3: EAF from study #2; ...
#' @param condp_cut Threshold of conditional p-value to be considered as significant. No default specified. Usually recommend 5e-8.
#' @param index_snps User-defined index SNP(s), if any. Default is `NULL` which means mJAM-Forward will automatically select index variants.
#' @param within_pop_threshold Threshold of r2 with selected index SNP(s) within a single population. If a SNP's correlation with any selected index SNP is greater than this threshold in at least one population, it will be excluded from subsequent rounds of index SNP selection.
#' @param across_pop_threshold Threshold of r2 with selected index SNP(s) across all populations. If a SNP's correlation with any selected index SNP is greater than this threshold in all populations, it will be excluded from subsequent rounds of index SNP selection.
#' @param coverage The required coverage of credible sets. Default is 0.95.
#' @param Pr_Med_cut Cut off of mJAM posterior mediation probability (P(Med)) during credible set construction. Low P(Med) may indicate low correlation between the candidate SNP and the index SNP. Any candidate credible set SNPs with P(Med) < Pr_Med_cut will be not be considered for credible set. Default is 0.
#' @param filter_rare A logical variable indicating whether to filter rare SNPs before the analysis. Default is `FALSE.` If `TRUE`, then please specify `rare_freq`.
#' @param rare_freq A vector of frequencies between 0 and 0.5 to specify the minor allele frequency cut-off if you want to filter rare SNPs before the analysis. Please also set `filter_rare` to be TRUE. For example, if there are 3 populations, then rare_freq = c(0.01, 0, 0.01) means SNPs with MAF < 0.01 in pop 1 and MAF < 0.01 in pop 3 will be removed from analysis.
#'
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
#'    \item{index}{A table listing all the selected index SNP(s) (`SNP`), along with their log10(p-value) conditional on all SNP(s) above (`cond_log10p`), the log10(p-value) conditional on all other index SNP(s) (`final_log10p`), and the p-value threshold used in this analysis (`pcut`).}
#'    \item{cs}{A table recording various posterior probabilities of all SNPs being considered for credible set SNPs. }
#'    \item{mJAM_marg_est}{A table with the marginal effect estimates and standard errors of all SNPs under the mJAM model.}
#'    \item{QC_marg_est}{The complete table of marginal effect estimates using fixed-effect model and mJAM model. For QC purpose only.}
#' }
#'


# N_GWAS = c(5000, 5000, 5000)
# X_ref = list(RefDosage_P1,RefDosage_P2,RefDosage_P3)
# Marg_Result = MargBeta
# EAF_Result = EAF
# condp_cut = 0.05/50
# index_snps = c("rs2")
# within_pop_threshold = 0.50
# across_pop_threshold = 0.20
# coverage = 0.95
# Pr_Med_cut = 0
# filter_rare = FALSE
# rare_freq = NULL

mJAM_Forward <- function(N_GWAS, X_ref,
                         Marg_Result, EAF_Result,
                         condp_cut = NULL,
                         index_snps = NULL,
                         within_pop_threshold = 0.50,
                         across_pop_threshold = 0.20,
                         coverage = 0.95,
                         Pr_Med_cut = 0,
                         filter_rare = FALSE,
                         rare_freq = NULL){
  ## Set parameters
  N_SNP <- numSNPs <- numSNPs_wo_rare <- nrow(Marg_Result)
  if(is.null(condp_cut)){condp_cut <- 0.05/N_SNP}

  ## Check index_snps is in marker names
  if(!is.null(index_snps) && sum(!(index_snps %in% Marg_Result$SNP))>0){
    not_in_index_snps = index_snps[!(index_snps %in% Marg_Result$SNP)]
    stop(paste0(paste0(not_in_index_snps, collapse = ","), "not found in Marg_Result."))
  }

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

  ## --- Calculate the LD matrix of X_ref
  Dosage_cor <- vector("list",length(X_ref))
  for(i in 1:length(X_ref)){
    if (nrow(missing_tbl) > 0 && i%in%missing_tbl$missing_ethnic_idx){
      ## --- Get missing SNP index
      temp_missing_snp_idx <- filter(missing_tbl, missing_ethnic_idx == i) %>% pull(missing_snp_idx)
      if(length(temp_missing_snp_idx)<ncol(X_ref[[i]])-1){
        ## --- Get Dosage_cor with complete SNP
        temp_Dosage_cor <- cor(X_ref[[i]][,-temp_missing_snp_idx])^2
        ## --- Fill in missing SNPs with zeros
        Dosage_cor[[i]] <- diag(1, nrow = numSNPs_wo_rare, ncol = numSNPs_wo_rare)
        Dosage_cor[[i]][-temp_missing_snp_idx, -temp_missing_snp_idx] <- temp_Dosage_cor
      }else{
        Dosage_cor[[i]] <- diag(1, nrow = numSNPs_wo_rare, ncol = numSNPs_wo_rare)
      }
    }else{
      Dosage_cor[[i]] <- cor(X_ref[[i]])^2
    }
  }

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

  ## Identify rare variants
  rare_pct_sumstats <- rowMeans(abs(EAF_Result[,2:ncol(EAF_Result)]-0.5)>0.48|is.na(EAF_Result[,2:ncol(EAF_Result)]))
  rare_SNPs_sumstats <- Marg_Result[which(rare_pct_sumstats>=0.5),1]
  reference_EAF = matrix(NA, nrow = numSNPs_wo_rare, ncol = length(X_ref))
  for(k in 1:length(X_ref)){
    reference_EAF[,k] <- colMeans(X_ref[[k]],na.rm=T)/2
  }
  rare_pct_dosage <- rowMeans(abs(reference_EAF-0.5)>0.48|is.na(reference_EAF))
  rare_SNPs_dosage <- Marg_Result[which(rare_pct_dosage>=0.5),1]
  rare_SNPs <- union(rare_SNPs_sumstats, rare_SNPs_dosage)
  if(length(rare_SNPs)==0){rare_SNPs <- NULL}

  ## Run Forward selection
  iter_count <- 0
  selected_ids <- NULL
  prev_index_snp <- NULL
  pruned_snps <- NULL
  condp_list <- NULL
  curr_min_condp <- 0
  prev_min_condp <- 0
  all_CS <- NULL
  subset_EUR <- has_dosage_SNP <- Marg_Result$SNP
  GItGI_curr <- GItGI
  GIty_curr <- GIty
  yty_curr <- yty

  while(iter_count >= 0){
    ## step 1: select top SNPs in the remaining list
    ## selected_id should be the ID in subset_EUR
    if(length(unique(pruned_snps))==numSNPs_wo_rare){break}
    if(iter_count == 0){Input_id = NULL}else{Input_id = match(selected_ids, subset_EUR)}

    ## get the id of rare SNPs in remaining SNPs
    if(sum(rare_SNPs %in% subset_EUR)>0){
      rare_id = match(rare_SNPs, subset_EUR)
      rare_id = rare_id[!is.na(rare_id)]
    }else{
      rare_id = NULL
    }

    ## get the conditional p-values of all remaining SNPs
    newFS_RES <- mJAM_get_condp(GItGI = GItGI_curr, GIty = GIty_curr, yty = yty_curr,
                                yty_med = yty_med, N_GWAS = N_GWAS, selected_id = Input_id,
                                use_median_yty_ethnic = NULL, rare_id = rare_id)

    ## output mJAM marginal p and meta marginal p
    if(iter_count == 0){
      marginal_est <- tibble(SNP = subset_EUR,
                             mJAM_effect = signif(newFS_RES$effect_est, digits = 3),
                             feMeta_effect = feMeta_mean,
                             mJAM_se = signif(newFS_RES$se_est, digits = 3),
                             feMeta_se = feMeta_se) %>%
        mutate(mJAM_logp = Rmpfr::pnorm(abs(mJAM_effect/mJAM_se),lower.tail = F,log.p = T)+log(2),
               feMeta_logp = Rmpfr::pnorm(abs(feMeta_effect/feMeta_se),lower.tail = F,log.p = T)+log(2)) %>%
        mutate(mJAM_log10p = signif(mJAM_logp/log(10), 4),
               feMeta_log10p = signif(feMeta_logp/log(10),4))
    }


    ## determine the index SNP of current round
    if(is.null(index_snps)){
      if((newFS_RES$condp_min>log(condp_cut))) {break}
      curr_index_snp <- subset_EUR[newFS_RES$which_condp_min]
    }else{
      if(iter_count >= length(index_snps)) {break}
      curr_index_snp <- index_snps[iter_count+1]
    }

    selected_ids <- c(selected_ids, curr_index_snp)
    condp_list <- c(condp_list, newFS_RES$condp[match(curr_index_snp, subset_EUR)])
    message(paste("No.", iter_count+1,"selected index SNP:", curr_index_snp))

    ## step 2: construct credible sets based on the selected SNP
    curr_CS <- mJAM_build_CS(X_id = curr_index_snp,
                             prev_X_list = prev_index_snp,
                             All_id = subset_EUR,
                             PrCS_weights = "Pr(M_C)",
                             coverage = coverage,
                             GItGI_curr = GItGI_curr, GIty_curr = GIty_curr,
                             yty_curr = yty_curr, yty_med = yty_med,
                             N_GWAS = N_GWAS, rare_SNPs = rare_SNPs,
                             Pr_Med_cut = Pr_Med_cut)

    ## step 3: prune out CS SNPs and SNPs in LD with index SNP; subset input statistics
    curr_CS_snp <- filter(curr_CS, CS_in == T) %>% pull(CS_SNP)
    pruned_snps <- c(pruned_snps, curr_index_snp)
    curr_LD_snp <- mJAM_LDpruning(target = match(curr_index_snp,has_dosage_SNP),
                                  testing = match(has_dosage_SNP[-match(pruned_snps,has_dosage_SNP)],has_dosage_SNP),
                                  R = Dosage_cor,
                                  within_thre = within_pop_threshold, across_thre = across_pop_threshold)
    pruned_snps <- c(pruned_snps,has_dosage_SNP[c(curr_LD_snp$remove_within, curr_LD_snp$remove_across)])

    all_CS <- rbind(all_CS, curr_CS)

    ## step 4: update input statistics
    if(length(pruned_snps)>0){
      subset_EUR <- has_dosage_SNP[-match(pruned_snps, has_dosage_SNP)]
    }
    subset_EUR <- union(subset_EUR, selected_ids)
    subset_EUR <- has_dosage_SNP[has_dosage_SNP%in%subset_EUR]
    message(paste(length(unique(pruned_snps))-length(unique(selected_ids)), "SNPs got pruned.", length(subset_EUR), "SNPs left."))

    if(length(subset_EUR)>1){
      for(e in 1:numEthnic){
        GItGI_curr[[e]] <- GItGI[[e]][subset_EUR, subset_EUR]
        GIty_curr[[e]] <- GIty[[e]][subset_EUR]
        yty_curr[[e]] <- yty[[e]][subset_EUR]
        colnames(GItGI_curr[[e]]) <- rownames(GItGI_curr[[e]]) <- subset_EUR
        names(GIty_curr[[e]]) <- names(yty_curr[[e]]) <- subset_EUR
      }
      ## step 5: continue to the next iteration
      iter_count <- iter_count+1
      prev_index_snp <- c(prev_index_snp, curr_index_snp)
      message("Continue to next round of index SNP selection.")
    }else{
      ## step 5: continue to the next iteration
      iter_count <- iter_count+1
      prev_index_snp <- c(prev_index_snp, curr_index_snp)
      break
    }

  }
  message("Analysis ended.")

  if(!is.null(selected_ids)){
    final_condp_selected <- mJAM_get_condp_selected(GItGI = GItGI, GIty = GIty, yty = yty,
                                                    N_GWAS = N_GWAS,yty_med = yty_med,
                                                    selected_id = selected_ids,
                                                    rare_SNPs = rare_SNPs)

    if(length(selected_ids)>1){
      finalp = final_condp_selected$condp[,1]
    }else{
      finalp =  final_condp_selected$condp
    }

    MULTI_index <- tibble(SNP = selected_ids,
                          b_joint = final_condp_selected$b_joint,
                          b_joint_var = ifelse(length(selected_ids)==1, final_condp_selected$b_joint_var,
                                               diag(final_condp_selected$b_joint_var)),
                          cond_log10p = condp_list/log(10),
                          final_log10p = finalp/log(10),
                          pcut = condp_cut)

    ## simplify CS output table
    all_CS <- all_CS %>%
      dplyr::select(c(index_SNP, CS_SNP, Post_Model_Prob_Ratio2, Post_Med_Prob2, SD_Post_CS_Prob, CumSum_Porb, CS_in)) %>%
      rename(Post_Model_Prob = Post_Model_Prob_Ratio2,
             Post_Med_Prob = Post_Med_Prob2,
             CumSum_Prob = CumSum_Porb)

  }else{
    message("No index SNP selected in this region.")
    MULTI_index = NULL
    all_CS = NULL
  }

  return(list(index = MULTI_index,
              cs = all_CS,
              mJAM_marg_est = marginal_est[,c(1,2,4,8)],
              QC_marg_est = marginal_est))
}


