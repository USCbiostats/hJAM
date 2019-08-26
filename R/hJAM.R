#' Fit hJAM
#' @description The hJAM function is to get the results from the hJAM model using input data
#'
#' @param betas.Gy The betas in the paper: the marginal effects of SNPs on the phenotype (Gy)
#' @param N.Gy The sample size of Gy
#' @param Gl The reference panel (Gl), such as 1000 Genome
#' @param Z The Z matrix in the paper: the marginal/conditional effects of SNPs on the exposures (Gx)
#' @param a_sigma The scale parameter of the prior, default = 1
#' @param b_sigma The scale parameter of the prior, default = 9
#' @param trait.variance The variance of the outcome, default = 1
#' @param ridgeTerm ridgeTerm = TRUE when the matrix L is singular. Matrix L is obtained from the cholesky decomposition of G0'G0. Default as FALSE.
#' @author Lai Jiang
#'
#' @examples
#' data(reference_data)
#' data(betas.Gy)
#' data(conditional_Z)
#' hJAM_sparse(betas.Gy = betas.Gy, Gl = Gl, N.Gy = 5000, Z = conditional_Z, ridgeTerm = TRUE)

hJAM_sparse = function(betas.Gy, N.Gy, Gl, Z, a_sigma = 1, b_sigma = 9, trait.variance = 1, ridgeTerm = FALSE) {

  # Check the dimension of betas.Gy, Gl and Z
  dim_betas = length(betas.Gy)
  dim_Gl = ncol(Gl)
  dim_Z = nrow(Z)

  if(dim_betas == dim_Gl & dim_betas == dim_Z){

    # The sample size in Gy
    N = N.Gy

    # Obtain the JAM variables: zL and L
    p = apply(Gl, 2, mean)/2
    n0 = N*(1-p)^2
    n1 = N*2*p*(1-p)
    n2 = N*p^2
    y0 = -(n1*betas.Gy+2*n2*betas.Gy)/(n0+n1+n2)
    y1 = y0+betas.Gy
    y2 = y0+2*betas.Gy
    z = n1*y1 + 2*n2*y2

    ## Compute G0'G0
    G0 = scale(Gl, center=T, scale=F)
    G0_t_G0 = t(G0)%*%G0

    ## Modify G0'G0 if the sample sizes of Gl and Gx are different
    Dm = 2*p*(1-p)*N
    D_sqrt = diag(sqrt(Dm))
    Dw_sqrt_inv = diag(1/sqrt(diag(G0_t_G0)))
    G0_t_G0.scaled = D_sqrt %*% Dw_sqrt_inv  %*% t(G0) %*% G0 %*% Dw_sqrt_inv %*% D_sqrt

    ## Add a ridge term in case G0'G0 is singular
    ridgeValue = ifelse(ridgeTerm, min(1,min(diag(G0_t_G0.scaled)*.001)), 0)
    G0_t_G0.ridge = G0_t_G0.scaled + ridgeValue*diag(length(betas.Gy))

    # For MR-JAM
    xTx = t(Z)%*% G0_t_G0.ridge %*% Z
    z = t(Z) %*% z
    g = N.Gy

    betas.XY = (solve(xTx) %*% z) *g/(1+g)

    # Get S2 for gibbs sampler
    # Old: s2 = t(zL-L%*%betas)%*%(zL-L%*%betas)
    # New IPD form (Bottolo 2010): y'y - 2X'yb +b'(X'X)b
    # y'y = trait.variance*(n.people-1)
    # X'y = z (NB not zL)
    s2 = trait.variance*(N.Gy-1) - 2*t(z)%*%betas.XY + (t(betas.XY)%*%xTx%*%betas.XY)

    # Estimate sigma2 (for use in Gibbs) - b/(a-1)
    # Old: pN = ncol(L); est_sigma2 = (b_sigma+s2/2+(t(betas)%*%t(L)%*%L%*%betas)/(2*(g+1)))/(a_sigma+pN/2-1)
    # n.people instead of pN.
    # Estimate looks MUCH better
    est_sigma2 = (b_sigma+s2/2+(t(betas.XY)%*%xTx%*%betas.XY)/(2*(g+1)))/(a_sigma+N.Gy/2-1)

    # Beta SEs in eq 13
    se_beta_temp = diag(solve(xTx))
    se.XY = sapply(se_beta_temp, function(x) sqrt(g*est_sigma2/(1+g)*x))

    pvalues.XY = 2*(1 - pnorm(abs(betas.XY/se.XY)))

    return(list(
      Estimate = betas.XY,
      StdErr = se.XY,
      Pvalue = pvalues.XY)
    )
  }else{
    cat("ERROR: The number of SNPs in betas.Gy, Z matrix and the reference panel (Gl) are different.")
  }
}
