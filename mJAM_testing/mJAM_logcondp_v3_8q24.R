# act_median <- 80000000
# yty <- c(.97, .93, .91, 6.5)*act_median
# med <- median(yty)
# delta <- abs(med-yty)
# w <- med/(med+delta)
# yty_tilde <- w*yty + (1-w)*med


mJAM_get_condp_selected = function(GItGI, GIty, yty,yty_med,N_GWAS, selected_id, 
                                   g = sum(N_GWAS), rare_SNPs = NULL){
  
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
    b_joint <- (g/(1+g))*susie_in_Xty/susie_in_XtX
    resid_var <- max(0,(susie_in_yty - susie_in_Xty*b_joint)/(sum_GWAS - length(selected_id)))
    b_joint_var <- (g/(1+g))*resid_var/susie_in_XtX
    sqrt(b_joint_var)
    condp_selected <- ifelse(resid_var==0, NA, Rmpfr::pnorm(abs(b_joint/sqrt(b_joint_var)), lower.tail = F, log.p = T)+ log(2))
  }else{
    # b_joint <- solve(susie_in_XtX) %*% susie_in_Xty
    # resid_var <- max(0,(susie_in_yty - t(susie_in_Xty)%*%b_joint)/(sum_GWAS - length(selected_id)))
    # b_joint_var <- resid_var*solve(susie_in_XtX)
    # condp_selected <- 2*pnorm(abs(b_joint/sqrt(diag(b_joint_var))), lower.tail = F)
    b_joint <- (g/(1+g))*solve(susie_in_XtX) %*% susie_in_Xty
    resid_var <- matrix(0,nrow = length(testing_ids), ncol = length(testing_ids))
    for(id in 1:length(testing_ids)){
      resid_var[id,id] <- max(0,(susie_in_yty[id] - t(susie_in_Xty)%*%b_joint)/(sum_GWAS[id] - length(testing_ids)))
    }
    b_joint_var <- (g/(1+g))*resid_var*solve(susie_in_XtX)
    sqrt_diag_b_joint_var <- ifelse(diag(b_joint_var)==0, NA, sqrt(diag(b_joint_var)))
    condp_selected <- Rmpfr::pnorm(abs(b_joint/sqrt_diag_b_joint_var), lower.tail = F, log.p = T)+ log(2)
  }
  
  return(list(b_joint = b_joint, b_joint_var = b_joint_var, cond_logp = condp_selected))
}


