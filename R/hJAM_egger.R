#hJAM_egger
#' Fit hJAM with Egger regression
#' @description The hJAM_egger function is to get the results from the hJAM model with Egger regression. It is for detecting potential pleiotropy
#'
#' @param betas.Gy The betas in the paper: the marginal effects of SNPs on the phenotype (Gy)
#' @param N.Gy The sample size of Gy
#' @param Gl The reference panel (Gl), such as 1000 Genome
#' @param A The A matrix in the paper: the marginal/conditional effects of SNPs on the exposures (Gx)
#' @param ridgeTerm ridgeTerm = TRUE when the matrix L is singular. Matrix L is obtained from the cholesky decomposition of G0'G0. Default as FALSE.
#' @author Lai Jiang
#'
#' @export
#'
#' @references
#'
#' Lai Jiang, Shujing Xu, Nicholas Mancuso, Paul J. Newcombe, David V. Conti (2020).
#' A Hierarchical Approach Using Marginal Summary Statistics for Multiple Intermediates
#' in a Mendelian Randomization or Transcriptome Analysis. \emph{bioRxiv}
#' \url{https://doi.org/10.1101/2020.02.03.924241}.
#'
#' @examples
#' data(Gl)
#' data(betas.Gy)
#' data(conditional_A)
#' hJAM_egger(betas.Gy = betas.Gy, Gl = Gl, N.Gy = 459324, A = conditional_A, ridgeTerm = TRUE)

hJAM_egger = function(betas.Gy, N.Gy, Gl, A, ridgeTerm = FALSE) {

  # Check the dimension of betas.Gy, Gl and A
  dim_betas = length(betas.Gy)
  dim_Gl = ncol(Gl)
  dim_A = nrow(A)

  if(dim_betas == dim_Gl & dim_betas == dim_A){

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
    X = cbind(rep(1, nrow(L)), L%*%A)
    betas.XY = summary(lm(zL ~ 0 + X))$coef[-1,1]
    se.XY = summary(lm(zL ~ 0 + X))$coef[-1,2]
    pvalues.XY = summary(lm(zL ~ 0 + X))$coef[-1,4]

    NaN_row = is.na(confint((lm(zL ~ 0 + X)))[, 1])
    lower.ci.all = confint((lm(zL ~ 0 + X)))[!NaN_row, 1]
    upper.ci.all = confint((lm(zL ~ 0 + X)))[!NaN_row, 2]

    lower.ci = lower.ci.all[-1]
    upper.ci = upper.ci.all[-1]

    betas.int = summary(lm(zL ~ 0 + X))$coef[1,1]
    se.int = summary(lm(zL ~ 0 + X))$coef[1,2]
    pvalues.int = summary(lm(zL ~ 0 + X))$coef[1,4]

    lower.ci.int = lower.ci.all[1]
    upper.ci.int = upper.ci.all[1]

    out <- list(
      Exposure = colnames(A),
      numSNP = nrow(X),
      Estimate = betas.XY,
      StdErr = se.XY,
      Pvalue = pvalues.XY,
      Lower.CI = lower.ci,
      Upper.CI = upper.ci,
      Est.Int = betas.int,
      StdErr.Int = se.int,
      Lower.CI.Int = lower.ci.int,
      Upper.CI.Int = upper.ci.int,
      Pvalue.Int = pvalues.int)
    class(out) <- "hJAM_egger"
    return(out)
  }else{
    cat("ERROR: The number of SNPs in betas.Gy, A matrix and the reference panel (Gl) are different.")
  }
}
