
#' Get Pr(Model) based on BF-type model probability
#'
#' @description Also apply weighting to get robust estimates of yty
#'
#' @param GItGI A list of transformed statistics from `get_XtX()` for each study.
#' @param GIty A list of transformed statistics from `get_z()` for each study.
#' @param yty A list of transformed statistics from `get_yty()` for each study.
#' @param yty_med A numeric vector of median yty across all SNPs within each study.
#' @param N_GWAS A numeric vector of GWAS sample size for each study.
#' @param C_id An ingeter vector of IDs for the SNPs to be tested.
#' @param prev_X_list A numeric vector of the ID(s) of previously selected index SNP(s).
#' @param g The pre-specified `g` in `g`-prior formulation.
#' @param rare_SNPs A numeric vector of ID(s) for rare SNP(s) which we do not apply weighting. Instead, we use the individual estimate of yty for these SNPs for robustness.
#'
#'
#' @author Jiayi Shen
#'
#' @importFrom utils capture.output
#' @importFrom Rmpfr mpfr
#'
#' @returns
#'
#' \describe{
#'    \item{post_prob}{Posterior Pr(Model) for each SNPs in `C_id`.}
#'    \item{R2_est}{R2 estimates of every one-SNP model (one for each SNPs in `C_id`).}
#'    \item{n_miss}{An integer vector of how many studies have missing values for each SNP.}
#' }
#'

mJAM_get_PrM <- function(GItGI,GIty,yty,yty_med,N_GWAS, C_id,prev_X_list = NULL, g = NULL,rare_SNPs = NULL){

  numEthnic <- length(GItGI)
  if(is.null(g)){g <- sum(N_GWAS)}

  ## --- Check whether c(C_id,prev_X_list) has missingness in any studies
  missing_ethnic_idx <- NULL
  for(i in 1:numEthnic){
    temp_GItGI_sum <- sum(GItGI[[i]][,c(C_id,prev_X_list)])
    if(temp_GItGI_sum == 0){
      missing_ethnic_idx <- c(missing_ethnic_idx, i)
    }
  }

  sum_GWAS = sum(N_GWAS) - sum(N_GWAS[missing_ethnic_idx])

  if(is.null(prev_X_list)){
    ## --- Calculate model-specific X'X and y'y
    GItGI_C <- GIty_C <- yty_C <- delta <- w <- vector("list", numEthnic)
    for(i in 1:numEthnic){
      GItGI_C[[i]] <- GItGI[[i]][c(C_id,prev_X_list),c(C_id,prev_X_list)]
      GIty_C[[i]] <- GIty[[i]][c(C_id,prev_X_list)]
      ## robust yty estimation
      ## robust yty estimation (v4.1b)
      delta[[i]] <- abs(yty[[i]][C_id] - yty_med[[i]])
      if(C_id %in% rare_SNPs){
        w[[i]] <- 1
      }else{
        w[[i]] <- yty_med[[i]]/(yty_med[[i]]+delta[[i]])
      }
      yty_C[[i]] <- w[[i]]*yty[[i]][C_id] + (1-w[[i]])*yty_med[[i]]
      # yty_C[[i]] <- yty_med[[i]] # (v4.1a)
      # yty_C[[i]] <-yty[[i]][C_id] # (v4.1c)
    }
    susie_in_XtX <- Reduce("+", GItGI_C)
    susie_in_Xty <- Reduce("+", GIty_C)
    susie_in_yty <- Reduce("+", yty_C)

    ## --- Calculate multi-ethnic R2 and Pr(Model)
    b_joint <- Rmpfr::mpfr(susie_in_Xty / susie_in_XtX,50)
    R2_est <-  b_joint*susie_in_XtX*b_joint/susie_in_yty
    R2 <- ifelse(R2_est > 1, 1, R2_est)
    post_prob <- 1/(1+g*(1-R2_est))^(0.5*(sum(N_GWAS)-1))

  }else{
    ## --- Calculate model-specific X'X and y'y
    GItGI_L <- GIty_L <- yty_L <- GItGI_S <- GIty_S <- vector("list", numEthnic)
    for(i in 1:numEthnic){
      GItGI_L[[i]] <- GItGI[[i]][c(C_id,prev_X_list),c(C_id,prev_X_list)]
      GIty_L[[i]] <- GIty[[i]][c(C_id,prev_X_list)]
      GItGI_S[[i]] <- GItGI[[i]][prev_X_list,prev_X_list]
      GIty_S[[i]] <- GIty[[i]][prev_X_list]
      if(i %in% missing_ethnic_idx){
        yty_L[[i]] <- 0
      }else{
        yty_L[[i]] <- yty_med[[i]]
      }
    }
    susie_in_XtX_L <- Reduce("+", GItGI_L)
    susie_in_Xty_L <- Reduce("+", GIty_L)
    susie_in_XtX_S <- Reduce("+", GItGI_S)
    susie_in_Xty_S <- Reduce("+", GIty_S)
    susie_in_yty <- Reduce("+", yty_L)

    ## SSR of larger model
    D_L <- diag(diag(susie_in_XtX_L))
    b_marg_L <- solve(D_L) %*% susie_in_Xty_L
    b_joint_L <- solve(susie_in_XtX_L) %*% D_L %*% b_marg_L
    SSR_L <- as.numeric(t(b_joint_L) %*% D_L %*% b_marg_L)
    ## SSR of smaller model
    if(is.null(dim(susie_in_XtX_S))){
      D_S <- susie_in_XtX_S
      b_marg_S <- susie_in_Xty_S/D_S
      b_joint_S <- (1/susie_in_XtX_S)*D_S*b_marg_S
      SSR_S <- as.numeric(b_joint_S * D_S * b_marg_S)
    }else{
      D_S <- diag(diag(susie_in_XtX_S))
      b_marg_S <- solve(D_S) %*% susie_in_Xty_S
      b_joint_S <- solve(susie_in_XtX_S) %*% D_S %*% b_marg_S
      SSR_S <- as.numeric(t(b_joint_S) %*% D_S %*% b_marg_S)
    }
    ## get partial R2
    R2_est <-  Rmpfr::mpfr( (SSR_L-SSR_S)/(susie_in_yty-SSR_S), 50)
    R2 <- ifelse(R2_est > 1, 1, R2_est)
    post_prob <- 1/(1+g*(1-R2_est))^(0.5*(sum(N_GWAS)-1))
    post_prob <- post_prob[1]
  }

  # return(post_prob)
  return(list(post_prob=post_prob,R2_est = R2_est, n_miss = length(missing_ethnic_idx)))
}

#' Get Pr(Model) based on Wald-type model probability
#'
#' @description Also apply weighting to get robust estimates of yty
#'
#' @param GItGI A list of transformed statistics from `get_XtX()` for each study.
#' @param GIty A list of transformed statistics from `get_z()` for each study.
#' @param yty A list of transformed statistics from `get_yty()` for each study.
#' @param yty_med A numeric vector of median yty across all SNPs within each study.
#' @param N_GWAS A numeric vector of GWAS sample size for each study.
#' @param C_id An ingeter vector of IDs for the SNPs to be tested.
#' @param prev_X_list A numeric vector of the ID(s) of previously selected index SNP(s).
#' @param g The pre-specified `g` in `g`-prior formulation.
#' @param rare_SNPs A numeric vector of ID(s) for rare SNP(s) which we do not apply weighting. Instead, we use the individual estimate of yty for these SNPs for robustness.
#'
#'
#' @author Jiayi Shen
#'
#' @importFrom stats pnorm
#'
#' @returns A numeric vector of posterior Pr(Model) for each SNPs in `C_id`.
#'
#'

mJAM_get_PrM_Wald <- function(GItGI, GIty, yty, yty_med, N_GWAS, C_id,prev_X_list = NULL, g = NULL, rare_SNPs = NULL){

  ## C -> Y | previous index SNPs

  ## Specify the prior information
  if(is.null(g)){g <- sum(N_GWAS)}
  numEthnic <- length(GItGI)

  ## Calculate model-specific X'X and y'y
  GItGI_C <- GIty_C <- yty_C <- w <- delta <- vector("list", numEthnic)
  for(i in 1:numEthnic){
    GItGI_C[[i]] <- GItGI[[i]][c(C_id,prev_X_list),c(C_id,prev_X_list)]
    GIty_C[[i]] <- GIty[[i]][c(C_id,prev_X_list)]
    ## robust yty estimation (v4.1)
    delta[[i]] <- abs(yty[[i]][C_id] - yty_med[[i]])
    if(C_id %in% rare_SNPs){
      w[[i]] <- 1
    }else{
      w[[i]] <- yty_med[[i]]/(yty_med[[i]]+delta[[i]])
    }
    yty_C[[i]] <- w[[i]]*yty[[i]][C_id] + (1-w[[i]])*yty_med[[i]]
  }

  multi_GItGI_C <- Reduce("+", GItGI_C)
  multi_GIty_C <- Reduce("+", GIty_C)
  multi_yty_C <- Reduce("+", yty_C)

  C_Y_s2 <- multi_yty_C - t(multi_GIty_C) %*% solve(multi_GItGI_C) %*% multi_GIty_C
  C_Y_bhat <- solve(multi_GItGI_C) %*% multi_GIty_C
  C_Y_post_var_scalar <- g*(C_Y_s2- t(C_Y_bhat) %*% multi_GItGI_C %*% (-C_Y_bhat)/(g+1))/(sum(N_GWAS)*(g+1))
  C_Y_post_var <- as.numeric(C_Y_post_var_scalar)*solve(multi_GItGI_C)
  C_Y_post_mean <- g/(g+1)*C_Y_bhat

  post_prob <- 2*(0.5-pnorm(abs(C_Y_post_mean[1]/sqrt(C_Y_post_var[1,1])), lower.tail = F))

  return(post_prob)
}

#' Get Pr(Mediation) based on causal mediation models
#'
#' @description Also apply weighting to get robust estimates of yty
#'
#' @param GItGI A list of transformed statistics from `get_XtX()` for each study.
#' @param GIty A list of transformed statistics from `get_z()` for each study.
#' @param yty A list of transformed statistics from `get_yty()` for each study.
#' @param yty_med A numeric vector of median yty across all SNPs within each study.
#' @param N_GWAS A numeric vector of GWAS sample size for each study.
#' @param g The pre-specified `g` in `g`-prior formulation.
#' @param C_id An ingeter vector of IDs for the SNPs to be tested.
#' @param X_id An integer specifying the ID of the index SNP.
#' @param prev_X_list A numeric vector of the ID(s) of previously selected index SNP(s).
#'
#'
#' @author Jiayi Shen
#'
#' @returns
#'
#' \describe{
#'    \item{Post_Med_Prob}{Posterior Pr(Mediation) for each SNPs in C_id.}
#'    \item{Med_Effect_Size}{Posterior mediation effect size for each SNPs in C_id.}
#'    \item{Med_var_CX}{Posterior variance of mediation effect in models with both C and X.}
#'    \item{Med_var_C}{Posterior variance of mediation effect in models with C only.}
#' }
#'

mJAM_get_PrMed <- function(GItGI, GIty,yty, yty_med, N_GWAS, g = NULL,C_id, X_id, prev_X_list){

  ## C -> Y | previous index SNPs

  ## Specify the prior information
  if(is.null(g)){g <- sum(N_GWAS)}
  numEthnic <- length(GIty)

  ## --- Check whether c(C_id,prev_X_list) has missingness in any studies
  missing_ethnic_idx <- NULL
  for(i in 1:numEthnic){
    temp_GItGI_sum <- sum(GItGI[[i]][,c(C_id,prev_X_list)])
    if(temp_GItGI_sum == 0){
      missing_ethnic_idx <- c(missing_ethnic_idx, i)
    }
  }

  sum_GWAS = sum(N_GWAS) - sum(N_GWAS[missing_ethnic_idx])

  ## Calculate model-specific X'X and y'y
  GItGI_C <- GIty_C <- yty_C <- vector("list", numEthnic)
  for(i in 1:numEthnic){
    GItGI_C[[i]] <- GItGI[[i]][c(C_id,prev_X_list),c(C_id,prev_X_list)]
    GIty_C[[i]] <- GIty[[i]][c(C_id,prev_X_list)]
    if(i %in% missing_ethnic_idx){
      yty_C[[i]] <- 0
    }else{
      yty_C[[i]] <- yty_med[[i]]
    }
  }

  multi_GItGI_C <- Reduce("+", GItGI_C)
  multi_GIty_C <- Reduce("+", GIty_C)
  multi_yty_C <- Reduce("+", yty_C)

  if((!is.null(dim(multi_GItGI_C))) && matrixcalc::is.singular.matrix(multi_GItGI_C)){
    ridgeValue <-  min(1,min(diag(multi_GItGI_C)*.001))
    multi_GItGI_C <-  multi_GItGI_C + ridgeValue*diag(dim(multi_GItGI_C)[2])
    message("multi_GItGI_C is singular. Ridge term added.")
  }
  ## s2 = (n-k-1)\hat{\sigma}^2 where \hat{\sigma}^2 is the estimated residual variance
  C_Y_s2 <- multi_yty_C - t(multi_GIty_C) %*% solve(multi_GItGI_C) %*% multi_GIty_C
  C_Y_bhat <- solve(multi_GItGI_C) %*% multi_GIty_C
  C_Y_post_var_scalar <- g*(C_Y_s2- t(C_Y_bhat) %*% multi_GItGI_C %*% (-C_Y_bhat)/(g+1))/(sum_GWAS*(g+1))
  C_Y_post_var <- as.numeric(C_Y_post_var_scalar)*solve(multi_GItGI_C)
  C_Y_post_mean <- g/(g+1)*C_Y_bhat


  ## C -> X -> Y | previous index SNPs
  ## --- Check whether c(X_id, C_id, prev_X_list) has missingness in any studies
  missing_ethnic_idx <- NULL
  for(i in 1:numEthnic){
    temp_GItGI_sum <- sum(GItGI[[i]][,c(X_id, C_id, prev_X_list)])
    if(temp_GItGI_sum == 0){
      missing_ethnic_idx <- c(missing_ethnic_idx, i)
    }
  }

  sum_GWAS = sum(N_GWAS) - sum(N_GWAS[missing_ethnic_idx])

  ## Specify the prior information
  GItGI_CX <- GIty_CX <- yty_CX <- vector("list", numEthnic)
  for(i in 1:numEthnic){
    GItGI_CX[[i]] <- GItGI[[i]][c(X_id, C_id, prev_X_list),c(X_id, C_id, prev_X_list)]
    GIty_CX[[i]] <- GIty[[i]][c(X_id, C_id, prev_X_list)]
    if(i %in% missing_ethnic_idx){
      yty_CX[[i]] <- 0
    }else{
      yty_CX[[i]] <- yty_med[[i]]
    }
  }

  multi_GItGI_CX <- Reduce("+", GItGI_CX)
  multi_GIty_CX <- Reduce("+", GIty_CX)
  multi_yty_CX <- Reduce("+", yty_CX)

  if((!is.null(dim(multi_GItGI_CX))) && matrixcalc::is.singular.matrix(multi_GItGI_CX)){
    ridgeValue <-  min(1,min(diag(multi_GItGI_CX)*.001))
    multi_GItGI_CX <-  multi_GItGI_CX + ridgeValue*diag(dim(multi_GItGI_CX)[2])
    message("multi_GItGI_CX is singular. Ridge term added.")
  }
  CX_Y_s2 <- multi_yty_CX - t(multi_GIty_CX) %*% solve(multi_GItGI_CX) %*% multi_GIty_CX
  CX_Y_bhat <- solve(multi_GItGI_CX) %*% multi_GIty_CX
  CX_Y_post_var_scalar <- g*(CX_Y_s2- t(CX_Y_bhat) %*% multi_GItGI_CX %*% (-CX_Y_bhat)/(g+1))/(sum_GWAS*(g+1))
  CX_Y_post_var <- as.numeric(CX_Y_post_var_scalar)*solve(multi_GItGI_CX)
  CX_Y_post_mean <- g/(g+1)*CX_Y_bhat

  med_test_stat <- (C_Y_post_mean[1] - CX_Y_post_mean[2])/sqrt(C_Y_post_var[1,1] + CX_Y_post_var[2,2])
  Post_Med_Prob <- 2*(0.5-pnorm(abs(med_test_stat), lower.tail = F))

  return(list(Post_Med_Prob = Post_Med_Prob,
              Med_Effect_Size = C_Y_post_mean[1] - CX_Y_post_mean[2],
              Med_var_CX = CX_Y_post_var[2,2],
              Med_var_C = C_Y_post_var[1,1]))
}

#' Construct mJAM credible set based for selected index SNP
#'
#'
#' @param X_id A character specifying the ID of the index SNP; should be found in `All_id`.
#' @param prev_X_list A list of character vector of the ID(s) of previously selected index SNP(s).
#' @param All_id A list of character vector of the ID(s) of all SNP(s) remaining in the analysis, including all previously selected SNP(s) and the current index SNP.
#' @param PrCS_weights An option to specify what weights to apply on Pr(Med). Default is "Pr(M_C)".
#' @param coverage A number between 0 and 1 specifying the “coverage” of the estimated confidence sets.
#' @param GItGI_curr A list of GItGI statistics at the current stage (after pruning out SNPs correlated with previously selected index SNPs).
#' @param GIty_curr A list of GIty estimates of all remaining SNPs at the current stage (after pruning out SNPs correlated with previously selected index SNPs).
#' @param yty_curr A list of yty estimates of all remaining SNPs at the current stage (after pruning out SNPs correlated with previously selected index SNPs).
#' @param yty_med A list of median yty across all SNPs.
#' @param N_GWAS A vector of sample sizes in all original GWAS studies.
#' @param rare_SNPs A numeric vector of ID(s) for rare SNP(s) which we do not apply weighting. Instead, we use the individual estimate of yty for these SNPs for robustness.
#' @param Pr_Med_cut The cutoff for Pr(Mediation); SNPs with Pr(Mediation) smaller than this cutoff will be assigned a Pr(CS) = 0 and thus not included in the credible set for the current index
#'
#'
#' @author Jiayi Shen
#'
#' @returns A table with the following columns:
#'
#' \describe{
#'    \item{CS_SNP}{SNP name.}
#'    \item{Post_Model_Prob}{The posterior Pr(Model) of this SNP on its absolute scale.}
#'    \item{Post_Model_Prob_Ratio}{The posterior Pr(Model) of this SNP divided by the posterior Pr(Model) of index SNP. It should be <= 1.}
#'    \item{Post_Model_Prob_Ratio2}{If `Post_Model_Prob_Ratio` is greater than 1, set `Post_Model_Prob_Ratio2` to 1. Otherwise, it's the same as `Post_Model_Prob_Ratio`.}
#'    \item{Med_Effect_Size}{The posterior mediation effect size.}
#'    \item{Post_Med_Prob}{The posterior Pr(Mediation) of this SNP.}
#'    \item{Post_Med_Prob2}{If `Post_Med_Prob` is less than `Pr_Med_cut`, set `Post_Med_Prob2` to 0. Otherwise, it's the same as `Post_Med_Prob`.}
#'    \item{SD_Post_CS_Prob}{Standardized Pr(CS) where Pr(CS) = Pr(Model)*Pr(Mediation)}
#'    \item{CumSum_Porb}{The cumulative `SD_Post_CS_Prob`. Note that the table is ordered by descending `SD_Post_CS_Prob`.}
#'    \item{EmpiricalCut}{The empirical coverage of this CS (should be >= requested `coverage`).}
#'    \item{CS_in}{A logical variable indicating whether this CS_SNP is included in this CS or not.}
#'    \item{index_SNP}{The name of the index SNP.}
#' }
#'
mJAM_build_CS <- function(X_id, prev_X_list = NULL,All_id, PrCS_weights = "Pr(M_C)", coverage = 0.95,
                          GItGI_curr, GIty_curr, yty_curr, yty_med, N_GWAS,rare_SNPs = NULL, Pr_Med_cut = 0.1){

  ###############################################################################################
  ## --- Initiate all variables
  Testing_id <- All_id[!All_id %in% c(X_id, prev_X_list)]
  Post_Med_Prob <-  rep(NA, length = length(All_id)); names(Post_Med_Prob) <- All_id
  Med_Effect_Size <- rep(NA, length(All_id)); names(Med_Effect_Size) <- All_id
  Post_Model_Prob <- rep(NA, length = length(All_id)); names(Post_Model_Prob) <- All_id
  Post_Model_Prob_Ratio <- rep(NA, length = length(All_id)); names(Post_Model_Prob_Ratio) <- All_id
  Post_CS_Prob <- rep(NA, length = length(All_id)); names(Post_CS_Prob) <- All_id

  ## --- Fill post probs for index SNP
  Post_Med_Prob[X_id] <- 1
  Med_Effect_Size[X_id] <- 0
  if(PrCS_weights == "Pr(M_C)"){
    temp_PrM_res <- mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr,
                                 yty = yty_curr, yty_med = yty_med,
                                 N_GWAS = N_GWAS, C_id = X_id,
                                 prev_X_list = prev_X_list, g = sum(N_GWAS),
                                 rare_SNPs = rare_SNPs)
    Post_Model_Prob[X_id] <- Post_CS_Prob[X_id] <- sum_Post_CS_Porb <- Index_Model_Prob <- temp_PrM_res$post_prob
    Post_Model_Prob_Ratio[X_id] <- 1
  }
  if(PrCS_weights == "Pr(Wald)"){
    Post_Model_Prob[X_id] <- Post_CS_Prob[X_id] <- sum_Post_CS_Porb <- Index_Model_Prob <-
      mJAM_get_PrM_Wald(GItGI = GItGI_curr, GIty = GIty_curr, yty = yty_curr, yty_med = yty_med,
                        N_GWAS = N_GWAS, C_id = X_id, rare_SNPs = rare_SNPs,
                        prev_X_list = prev_X_list, g = sum(N_GWAS))
    Post_Model_Prob_Ratio[X_id] <- 1
  }


  for(C_id in Testing_id){
    temp_med <- mJAM_get_PrMed(GItGI = GItGI_curr, GIty = GIty_curr,yty = yty_curr,
                               yty_med = yty_med,
                               N_GWAS  = N_GWAS, g = sum(N_GWAS),
                               C_id = C_id,
                               X_id = X_id,
                               prev_X_list = NULL)
    Post_Med_Prob[C_id] <- temp_med$Post_Med_Prob
    Med_Effect_Size[C_id] <- temp_med$Med_Effect_Size

    if(PrCS_weights == "Pr(M_C)"){
      temp_PrM_res <- mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr,
                                   yty = yty_curr, yty_med = yty_med,
                                   N_GWAS = N_GWAS, C_id = C_id,
                                   prev_X_list = prev_X_list, g = sum(N_GWAS),
                                   rare_SNPs = rare_SNPs)
      temp_r2 <- temp_PrM_res$post_prob
    }
    if(PrCS_weights == "Pr(Wald)"){
      temp_r2 <- mJAM_get_PrM_Wald(GItGI = GItGI_curr, GIty = GIty_curr, N_GWAS = N_GWAS,
                                   yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs,
                                   C_id = C_id, prev_X_list = prev_X_list, g = sum(N_GWAS))
    }
    Post_Model_Prob[C_id] <- temp_r2
    Post_Model_Prob_Ratio[C_id] <- as.double(temp_r2 / Index_Model_Prob)
    Post_CS_Prob[C_id] <- temp_r2*Post_Med_Prob[C_id]
    sum_Post_CS_Porb <- sum_Post_CS_Porb + temp_r2*Post_Med_Prob[C_id]
  }

  ###############################################################################################

  ## --- Get "standardized" posterior model probability
  Post_Model_Prob <- rep(NA, length = length(All_id)); names(Post_Model_Prob) <- All_id
  # Post_CS_Prob <- rep(NA, length = length(All_id)); names(Post_CS_Prob) <- All_id
  Std_CS_Prob <- rep(NA, length(All_id)); names(Std_CS_Prob) <- All_id

  ## --- Fill post probs for index SNP
  if(PrCS_weights == "Pr(Wald)"){
    Post_Med_Prob[X_id] <- 1
    Post_Model_Prob[X_id] <- Post_CS_Prob[X_id] <-
      mJAM_get_PrM_Wald(GItGI = GItGI_curr, GIty = GIty_curr,
                        yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs,
                        N_GWAS = N_GWAS, C_id = X_id, prev_X_list = prev_X_list,
                        g = sum(N_GWAS))
    Std_CS_Prob[X_id] <- Post_Model_Prob[X_id]/sum_Post_CS_Porb
  }
  if(PrCS_weights == "Pr(M_C)"){
    Post_Med_Prob[X_id] <- 1
    Post_Model_Prob[X_id] <- Post_CS_Prob[X_id] <-
      mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr,
                   yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs,
                   N_GWAS = N_GWAS, C_id = X_id, prev_X_list = prev_X_list,
                   g = sum(N_GWAS))$post_prob
    Std_CS_Prob[X_id] <-
      as.numeric(mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr,
                              yty = yty_curr, yty_med = yty_med,rare_SNPs = rare_SNPs,
                              N_GWAS = N_GWAS, C_id = X_id, prev_X_list = prev_X_list,
                              g = sum(N_GWAS))$post_prob/sum_Post_CS_Porb)
  }

  for(C_id in Testing_id){

    if(PrCS_weights == "Pr(Wald)"){
      temp_r2 <-  mJAM_get_PrM_Wald(GItGI = GItGI_curr, GIty = GIty_curr,
                                    yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs,
                                    N_GWAS = N_GWAS, C_id = C_id, prev_X_list = prev_X_list,
                                    g = sum(N_GWAS))
    }
    if(PrCS_weights == "Pr(M_C)"){
      temp_r2 <- mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr,
                              yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs,
                              N_GWAS = N_GWAS, C_id = C_id, prev_X_list = prev_X_list,
                              g = sum(N_GWAS))$post_prob
    }

    Post_Model_Prob[C_id] <- temp_r2
    Std_CS_Prob[C_id] <- as.numeric(temp_r2*Post_Med_Prob[C_id]/sum_Post_CS_Porb)
  }

  ###############################################################################################
  ## 2021.09.29 add in Post_Model_Prob, Post_Model_Prob_Ratio, Post_Med_Prob, Med_Effect_Size
  if(PrCS_weights == "Pr(M_C)"){
    Post_Model_Prob_str <- rep(NA, length(Post_Model_Prob))
    for(j in 1:length(Post_Model_Prob)){
      Post_Model_Prob_str[j] <- sub("\'mpfr1\' ", "", capture.output(Post_Model_Prob[[names(Post_Model_Prob[j])]]))
    }
  }
  if(PrCS_weights == "Pr(Wald)"){
    Post_Model_Prob_str <- Post_Model_Prob
  }

  Post_Model_Prob_Ratio2 <- Post_Model_Prob_Ratio
  Post_Model_Prob_Ratio2[which(Post_Model_Prob_Ratio>1)] <- 1
  Post_Med_Prob2 <- Post_Med_Prob
  Post_Med_Prob2[which(Post_Med_Prob<Pr_Med_cut)] <- 0
  SD_Post_CS_Prob2 <- (Post_Model_Prob_Ratio2*Post_Med_Prob2)/sum(Post_Model_Prob_Ratio2*Post_Med_Prob2, na.rm = T)

  Post_Med_Prob_df <-
    tibble(SNP_names = All_id,
           Post_Model_Prob = Post_Model_Prob_str,
           Post_Model_Prob_Ratio = Post_Model_Prob_Ratio,
           Post_Model_Prob_Ratio2 = Post_Model_Prob_Ratio2,
           Med_Effect_Size = Med_Effect_Size,
           Post_Med_Prob = Post_Med_Prob,
           Post_Med_Prob2 = Post_Med_Prob2,
           SD_Post_CS_Prob = SD_Post_CS_Prob2) %>%
    arrange(desc(SD_Post_CS_Prob)) %>%
    mutate(CumSum_Porb = cumsum(SD_Post_CS_Prob),
           EmpiricalCut = min(CumSum_Porb[CumSum_Porb >= coverage], na.rm = T),
           duplicated_one = duplicated(CumSum_Porb),
           CS_in = ((CumSum_Porb <= EmpiricalCut) & !duplicated_one))%>%
    dplyr::select(-duplicated_one)

  temp_CS_set <- Post_Med_Prob_df %>%
    rename(CS_SNP = SNP_names) %>%
    mutate(index_SNP = X_id)

  return(temp_CS_set)
}

