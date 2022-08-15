
#' Run mJAM with SuSiE
#'
#' @description fitting mJAM-SuSiE
#'
#' @param marginal.betas A list of marginal effect estimates from each study/population. The length of the list should be equal to the number of studies.
#' @param marginal.se A list of marginal standard error of `marginal.betas` from each study/population.
#' @param EAFs A list of effect allele frequencies from each study/population.
#' @param N_GWAS Sample sizes in all original GWAS studies.
#' @param X_ref A list of matrices with individual-level SNP dosage data in each study/population.
#' @param SNP_names A character vector of the names of the SNPs. Note that the order of SNP names should be the same as `EAFs`, `X_ref`, `marginal.betas` and `marginal.se`.
#' @param logORToBeta If `TRUE`, transform the marginal logOR estimates in `marginal.betas` to marginal linear effects.
#' @param p_cases (Required if `logORToBeta == TRUE`) The proportion of cases in each study/population, use `c()` to combine.
#' @param replace_missing If `TRUE`, then replace missing SNPs with zeros in the summary statistics of the missing study(s). Default is TRUE.
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

  # if("susieR" %in% rownames(installed.packages()) == FALSE) {install.packages("susieR")}
  # library(susieR); library(dplyr); library(tibble)
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
