
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
  n0 = N_outcome*(1-maf)
  n1 = N_outcome*maf

  y0 = -(n1*betas)/(n0+n1)
  y1 = y0+betas
  z = n1*y1

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
    Dj <- maf*(1-maf)*N_outcome
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
    Dj <- maf*(1-maf)*N_outcome
    D_sqrt <- sqrt(Dj)
    Dw_sqrt_inv <- 1/sqrt(G0_t_G0)
    G0_t_G0.scaled <- D_sqrt * Dw_sqrt_inv * G0_t_G0 * Dw_sqrt_inv * D_sqrt
  }

  return(G0_t_G0.scaled)
}


#' Get transformed statistics: yty
#'
#' @description
#' To calculate sufficient statistics based on summary statistics.
#' This yty estimate follows Yang et al. (2012) Nat Gen.
#' Marginal estimates from one SNP will produce one yty estimates.
#' Yang suggests taking the median across all SNPs to obtain a robust estimate.
#' Here we record all yty estimates and output both the median and the entire vector.
#'
#' @param maf A vector of minor allele frequencies
#' @param N_outcome Sample size in the GWAS where we obtained `betas`
#' @param betas A vector of marginal estimates of effect sizes (betas for continuous outcome; logOR for binary outcome)
#' @param betas.se A vector of the standard errors of marginal effect estimates (`betas`).

#'
#' @returns median of yty estimates across all SNPs; and a vector of all yty estimates

get_yty <- function(maf, N_outcome, betas, betas.se){

  ## Follows Yang's Nature paper (2012)
  Dj <- maf*(1-maf)*N_outcome
  Sj2 <- betas.se^2
  yTy.all <- Dj*Sj2*(N_outcome-1)+Dj*betas^2
  yTy.est <- median(Dj*Sj2*(N_outcome-1)+Dj*betas^2, na.rm = TRUE)
  # yTy.all[is.na(yTy.all)] <- median(Dj*Sj2*(N_outcome-1)+Dj*betas^2, na.rm = TRUE)
  yTy.all[is.na(yTy.all)] <-  0

  return(list(yTy.est = yTy.est, yTy.all = yTy.all))
}


#' Transform log odds ratios to linear effects
#'
#' @description
#' Adopted from R2BGLiMS::JAM_LogisticToLinearEffects. Reference: Benner 2015, FINEMAP
#'
#' @param log.ors A vector of log odds ratios
#' @param log.or.ses A vector of the standard errors of the log ORs
#' @param snp.genotype.sds A vector of standard deviations of genotypes (optional if `mafs` is provided)
#' @param mafs A vector of effective allele frequencies (optional if `snp.genotype.sds` is provided)
#' @param n Sample size in the GWAS where we obtained `log.ors`
#' @param p.cases A numeric value of the proportion of cases in the GWAS.
#'
#' @returns Transformed linear effect estimates, and transformed standards errors of linear effects.

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

