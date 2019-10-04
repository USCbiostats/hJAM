#hJAM_egger
#' Fit hJAM with Egger regression
#' @description The hJAM_egger function is to get the results from the hJAM model with Egger regression. It is for detecting potential pleiotropy
#'
#' @param betas.Gy The betas in the paper: the marginal effects of SNPs on the phenotype (Gy)
#' @param N.Gy The sample size of Gy
#' @param Gl The reference panel (Gl), such as 1000 Genome
#' @param Z The Z matrix in the paper: the marginal/conditional effects of SNPs on the exposures (Gx)
#' @param ridgeTerm ridgeTerm = TRUE when the matrix L is singular. Matrix L is obtained from the cholesky decomposition of G0'G0. Default as FALSE.
#' @author Lai Jiang
#'
#' @export
#' @importFrom stats lm pnorm
#' @examples
#' data(reference_data)
#' data(betas.Gy)
#' data(conditional_Z)
#' hJAM_egger(betas.Gy = betas.Gy, Gl = Gl, N.Gy = 5000, Z = conditional_Z, ridgeTerm = TRUE)

hJAM_egger = function(betas.Gy, N.Gy, Gl, Z, ridgeTerm = FALSE) {

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
    G0_t_G0.scaled = D_sqrt %*% Dw_sqrt_inv  %*% G0_t_G0 %*% Dw_sqrt_inv %*% D_sqrt

    ## Add a ridge term in case G0'G0 is singular
    ridgeValue = ifelse(ridgeTerm, min(1,min(diag(G0_t_G0.scaled)*.001)), 0)
    G0_t_G0.ridge = G0_t_G0.scaled + ridgeValue*diag(length(betas.Gy))

    # Perfrom Cholesky decompostion and construct zL
    L = chol(G0_t_G0.ridge)
    zL = solve(t(L))%*%z

    # Perform linear regression
    X = cbind(rep(1, nrow(L)), L%*%Z)
    betas.XY = summary(lm(zL ~ 0 + X))$coef[-1,1]
    se.XY = summary(lm(zL ~ 0 + X))$coef[-1,2]
    pvalues.XY = summary(lm(zL ~ 0 + X))$coef[-1,4]

    betas.int = summary(lm(zL ~ 0 + X))$coef[1,1]
    se.int = summary(lm(zL ~ 0 + X))$coef[1,2]
    pvalues.int = summary(lm(zL ~ 0 + X))$coef[1,4]

    return(list(
      Estimate = betas.XY,
      StdErr = se.XY,
      Pvalue = pvalues.XY,
      Est.Int = betas.int,
      StdErr.Int = se.int,
      Pvalue.Int = pvalues.int)
    )
  }else{
    cat("ERROR: The number of SNPs in betas.Gy, Z matrix and the reference panel (Gl) are different.")
  }
}
