
#' Get transformed statistics: z, or Xty
#'
#' @description To calculate sufficient statistics based on summary statistics
#'
#' @param maf A vector of minor allele frequencies
#' @param betas A vector of marginal estimates of effect sizes (betas for continuous outcome; logOR for binary outcome)
#' @param N_outcome Sample size in the GWAS where we obtained `betas`
#'
#' @returns a numeric vector of calculated z statistic

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


#' Get transformed statistics: XtX
#'
#' @description To calculate sufficient statistics based on summary statistics
#'
#' @param N_outcome Sample size in the GWAS where we obtained `betas`
#' @param Gl A matrix of reference dosage, columns are SNPs and rows are individuals.
#' @param maf A vector of minor allele frequencies
#'
#' @returns a variance covariance matrix of scaled Gl

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

