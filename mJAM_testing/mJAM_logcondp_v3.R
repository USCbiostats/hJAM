# act_median <- 80000000
# yty <- c(.97, .93, .91, 6.5)*act_median
# med <- median(yty)
# delta <- abs(med-yty)
# w <- med/(med+delta)
# yty_tilde <- w*yty + (1-w)*med


mJAM_get_condp_selected = function(GItGI, GIty, yty,yty_med,N_GWAS, selected_id, 
                                   use_median_yty_ethnic = NULL, rare_SNPs = NULL){
  
  ## --- Check whether selected_id has missingness in any studies
  numEthnic <- length(N_GWAS)
  testing_ids <- selected_id
  missing_ethnic_idx <- vector("list", length(testing_ids))
  sum_GWAS <- rep(sum(N_GWAS),length(testing_ids))
  for(id in 1:length(testing_ids)){
    for(i in 1:numEthnic){
      temp_GItGI_sum <- sum(GItGI[[i]][,testing_ids[id]])
      if(temp_GItGI_sum == 0){
        missing_ethnic_idx[[id]] <- c(missing_ethnic_idx[[id]], i)
      }
    }
    sum_GWAS[id] <- sum(N_GWAS) - sum(N_GWAS[missing_ethnic_idx[[id]]])
  }
  
  
  ## --- Get model-specific sufficient statistics
  GItGI_sub <- GIty_sub <- yty_sub <- yty_ind <- yty_median <-  vector("list", numEthnic)
  w <- delta <- vector("list", numEthnic)
  for(i in 1:numEthnic){
    GItGI_sub[[i]] <- GItGI[[i]][testing_ids,testing_ids]
    GIty_sub[[i]] <- GIty[[i]][testing_ids]
    yty_ind[[i]] <- yty[[i]][testing_ids]
    yty_median[[i]] <- rep(yty_med[[i]], length(testing_ids))
    delta[[i]] <- abs(yty_ind[[i]] - yty_median[[i]])
    w[[i]] <- yty_median[[i]]/(yty_median[[i]]+delta[[i]])
    w[[i]][which(testing_ids %in% rare_SNPs)] <- 1
    yty_sub[[i]] <- w[[i]]*yty_ind[[i]] + (1-w[[i]])*yty_median[[i]]
  }
  
  ## replace yty estimate with median in ethnic groups with large variability in yty estimates
  if(!is.null(use_median_yty_ethnic)){
    for(rep_id in use_median_yty_ethnic){
      yty_sub[[rep_id]] <- rep(median(yty[[rep_id]], na.rm = T), length(testing_ids))
    }
  }
  
  ## for missing SNPs, fill in yty with 0
  for(id in 1:length(testing_ids)){
    temp_missing_id <- missing_ethnic_idx[[id]]
    for(ethnic_id in temp_missing_id){
      yty_sub[[ethnic_id]][id] <- 0
    }
  }
  
  susie_in_XtX <- Reduce("+", GItGI_sub)
  susie_in_Xty <- Reduce("+", GIty_sub)
  susie_in_yty <- Reduce("+", yty_sub)
  
  ## --- Calculate multi-ethnic R-squared
  if(is.null(dim(susie_in_XtX))){
    b_joint <- susie_in_Xty / susie_in_XtX
    resid_var <- max(0,(susie_in_yty - susie_in_Xty*b_joint)/(sum_GWAS - length(selected_id)))
    b_joint_var <- resid_var/susie_in_XtX
    sqrt(b_joint_var)
    condp_selected <- ifelse(resid_var==0, NA, Rmpfr::pnorm(abs(b_joint/sqrt(b_joint_var)), lower.tail = F, log.p = T)+ log(2))
  }else{
    # b_joint <- solve(susie_in_XtX) %*% susie_in_Xty
    # resid_var <- max(0,(susie_in_yty - t(susie_in_Xty)%*%b_joint)/(sum_GWAS - length(selected_id)))
    # b_joint_var <- resid_var*solve(susie_in_XtX)
    # condp_selected <- 2*pnorm(abs(b_joint/sqrt(diag(b_joint_var))), lower.tail = F)
    b_joint <- solve(susie_in_XtX) %*% susie_in_Xty
    resid_var <- matrix(0,nrow = length(testing_ids), ncol = length(testing_ids))
    for(id in 1:length(testing_ids)){
      resid_var[id,id] <- max(0,(susie_in_yty[id] - t(susie_in_Xty)%*%b_joint)/(sum_GWAS[id] - length(testing_ids)))
    }
    b_joint_var <- resid_var*solve(susie_in_XtX)
    sqrt_diag_b_joint_var <- ifelse(diag(b_joint_var)==0, NA, sqrt(diag(b_joint_var)))
    condp_selected <- Rmpfr::pnorm(abs(b_joint/sqrt_diag_b_joint_var), lower.tail = F, log.p = T)+ log(2)
  }
  
  return(list(b_joint = b_joint, b_joint_var = b_joint_var, condp = condp_selected))
}


mJAM_get_condp = function(GItGI, GIty, yty, yty_med, N_GWAS, selected_id, use_median_yty_ethnic = NULL, rare_id = NULL){
  
  ## --- Specify which SNPs to include
  numSNPs <- length(GIty[[1]])
  numEthnic <- length(N_GWAS)
  condp <- rep(NA, numSNPs)
  effect_est <- rep(NA, numSNPs)
  se_est <- rep(NA, numSNPs)
  p_model <- length(selected_id)+1
  condp_mx <- matrix(NA,nrow = numSNPs, ncol = 3*p_model)
  All_SNPs <- 1:numSNPs
  in_id_list <- All_SNPs[!All_SNPs %in% selected_id]
  
  for (in_id in in_id_list){
    ## --- Check whether c(in_id, selected_id) has missingness in any studies
    testing_ids <- c(in_id, selected_id)
    missing_ethnic_idx <- vector("list", length(testing_ids))
    sum_GWAS <- rep(sum(N_GWAS),length(testing_ids))
    for(id in 1:length(testing_ids)){
      for(i in 1:numEthnic){
        temp_GItGI_sum <- sum(GItGI[[i]][,testing_ids[id]])
        if(temp_GItGI_sum == 0){
          missing_ethnic_idx[[id]] <- c(missing_ethnic_idx[[id]], i)
        }
      }
      sum_GWAS[id] <- sum(N_GWAS) - sum(N_GWAS[missing_ethnic_idx[[id]]])
    }
    
    ## --- Get model-specific sufficient statistics
    GItGI_sub <- GIty_sub <- yty_sub <- yty_ind <- yty_median <-  vector("list", numEthnic)
    w <- delta <- vector("list", numEthnic)
    for(i in 1:numEthnic){
      GItGI_sub[[i]] <- GItGI[[i]][testing_ids,testing_ids]
      GIty_sub[[i]] <- GIty[[i]][testing_ids]
      yty_ind[[i]] <- yty[[i]][testing_ids]
      yty_median[[i]] <- rep(yty_med[[i]], length(testing_ids))
      delta[[i]] <- abs(yty_ind[[i]] - yty_median[[i]])
      # if(i==4){
      #   w[[i]] <- yty_median[[i]]/(yty_median[[i]]+delta[[i]])
      # }else{
      #   w[[i]] <- 0
      # }
      # w[[i]] <- yty_median[[i]]/(yty_median[[i]]+delta[[i]])
      w[[i]] <- yty_median[[i]]/(yty_median[[i]]+delta[[i]])
      w[[i]][which(testing_ids %in% rare_id)] <- 1
      yty_sub[[i]] <- w[[i]]*yty_ind[[i]] + (1-w[[i]])*yty_median[[i]]
    }
    
    
    
    ## replace yty estimate with median in ethnic groups with large variability in yty estimates
    if(!is.null(use_median_yty_ethnic)){
      for(rep_id in use_median_yty_ethnic){
        yty_sub[[rep_id]] <- rep(median(yty[[rep_id]], na.rm = T), length(testing_ids))
      }
    }
    
    
    for(id in 1:length(testing_ids)){
      temp_missing_id <- missing_ethnic_idx[[id]]
      for(ethnic_id in temp_missing_id){
        yty_sub[[ethnic_id]][id] <- 0
      }
    }
    
    susie_in_XtX <- Reduce("+", GItGI_sub)
    susie_in_Xty <- Reduce("+", GIty_sub)
    susie_in_yty <- Reduce("+", yty_sub)
    
    ## --- Calculate multi-ethnic R-squared
    if(is.null(dim(susie_in_XtX))){
      b_joint <- susie_in_Xty / susie_in_XtX
      resid_var <- max(0,(susie_in_yty - susie_in_Xty*b_joint)/(sum_GWAS - length(testing_ids)))
      b_joint_var <- resid_var/susie_in_XtX
      condp[in_id] <- ifelse(resid_var==0, NA, Rmpfr::pnorm(abs(b_joint/sqrt(b_joint_var)), lower.tail = F, log.p = T)+log(2))
      effect_est[in_id] <- b_joint
      se_est[in_id] <- sqrt(b_joint_var)
      condp_mx[in_id,1] <- b_joint
      condp_mx[in_id,2] <- sqrt(b_joint_var)
      condp_mx[in_id,3] <- ifelse(resid_var==0, NA, Rmpfr::pnorm(abs(b_joint/sqrt(b_joint_var)), lower.tail = F, log.p = T)+log(2))
    }else{
      b_joint <- solve(susie_in_XtX) %*% susie_in_Xty
      resid_var <- matrix(0,nrow = length(testing_ids), ncol = length(testing_ids))
      for(id in 1:length(testing_ids)){
        resid_var[id,id] <- max(0,(susie_in_yty[id] - t(susie_in_Xty)%*%b_joint)/(sum_GWAS[id] - length(testing_ids)))
      }
      b_joint_var <- resid_var*solve(susie_in_XtX)
      condp[in_id] <- ifelse(resid_var[1,1]==0, NA, Rmpfr::pnorm(abs(b_joint[1]/sqrt(b_joint_var[1,1])), lower.tail = F, log.p = T)+log(2))
      effect_est[in_id] <- b_joint[1,1]
      se_est[in_id] <- sqrt(b_joint_var[1,1])
      D <- diag(diag(susie_in_XtX))
      b_marg <- solve(D) %*% susie_in_Xty
      condp_mx[in_id, 3*(0:length(selected_id))+1] <- as.vector(b_joint)
      condp_mx[in_id, 3*(0:length(selected_id))+2] <- sqrt(diag(b_joint_var))
      condp_mx[in_id, 3*(0:length(selected_id))+3] <- Rmpfr::pnorm(abs(as.vector(b_joint)/sqrt(diag(b_joint_var))),lower.tail = F, log.p = T)+log(2)
    }
  }
  
  # tibble(zscore = zscore, R2 = R2) %>%
  #   ggplot(aes(x = zscore, y = R2)) +
  #   geom_point(color="grey")
  # ggsave(paste0(Forward_Dir,"R2_diag_plot.png"), width = 8, height = 5)
  
  return(list(which_condp_min =  which.min(condp), 
              condp_min = condp[ which.min(condp)], 
              condp = condp, 
              effect_est = effect_est, se_est = se_est,
              condp_mx = condp_mx
  ))
}
