# MetaXcan function
metaXcan = function(betas, se.betaG, Gamma_G, Z, sigma_g, sigma_l){
  if(dim(sigma_g)[1] != 1){
    cat('dim(sigma_g)[1] != 1 --> More than one X in the model \n')
    break
  }else{
    beta_metaXcan = sum(Z*betas*(sigma_l^2))/(sigma_g^2) #beta to return
    vec.sigma_G = rep(sigma_g, times=length(sigma_l))
    z.score = sum(Z*(sigma_l/vec.sigma_G)*(betas/se.betaG))
    p = 2*pnorm(-abs(z.score))
    se_metaXcan = beta_metaXcan/z.score
    return(list(betas=beta_metaXcan, ses=se_metaXcan, p=p))
  }
}
