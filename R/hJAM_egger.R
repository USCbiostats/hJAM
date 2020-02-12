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
#' @return An object of the hJAM with egger regression results.
#' \describe{
#' \item{Exposure}{The intermediates, such as the modifiable risk factors in Mendelian Randomization and
#' gene expression in transcriptome analysis.}
#' \item{numSNP}{The number of SNPs that the user use in the instrument set.}
#' \item{Estimate}{The conditional estimates of the associations between intermediates and the outcome.}
#' \item{StdErr}{The standard error of the conditional estimates of the associations between intermediates
#' and the outcome.}
#' \item{Lower.CI}{The lower bound of the 95\% confidence interval of the estimates.}
#' \item{Upper.CI}{The upper bound of the 95\% confidence interval of the estimates.}
#' \item{Pvalue}{The p value of the estimates with a type-I error equals 0.05.}
#' \item{Est.Int}{The intercept of the regression of intermediates on the outcome.}
#' \item{StdErr.Int}{The standard error of the intercept of the regression of intermediates
#' on the outcome.}
#' \item{Lower.CI.Int}{The lower bound of the 95\% confidence interval of the intercept.}
#' \item{Upper.CI.Int}{The upper bound of the 95\% confidence interval of the intercept.}
#' \item{Pvalue.Int}{The p value of the intercept with a type-I error equals 0.05.}
#' }
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

#' @return An object of hJAM with egger regression results.

hJAM_egger = function(betas.Gy, N.Gy, Gl, A, ridgeTerm = FALSE) {

  # Check the dimension of betas.Gy, Gl and A
  dim_betas = length(betas.Gy)
  dim_Gl = ncol(Gl)
  dim_A = ifelse(is.null(dim(A)), length(A), nrow(A))

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
    glm.out = summary(glm(zL ~ 0 + X, family = gaussian()))
    betas.XY = glm.out$coef[-1,1]
    se.XY = glm.out$coef[-1,2]
    pvalues.XY = 2*pnorm(-abs(betas.XY/se.XY))

    lower.ci = betas.XY+qnorm(0.05)*se.XY
    upper.ci = betas.XY+qnorm(0.95)*se.XY

    betas.int = glm.out$coef[1,1]
    se.int = glm.out$coef[1,2]
    pvalues.int = 2*pnorm(-abs(betas.int/se.int))

    lower.ci.int = betas.int+qnorm(0.05)*se.int
    upper.ci.int = betas.int+qnorm(0.95)*se.int

    if(is.null(colnames(A))){
      if(!is.null(dim(A))){
        colnames_A = paste0("exposure ", 1:ncol(A))
      }else{
        colnames_A = "exposure"
      }
    }else{
      colnames_A = colnames(A)
    } # add column names of A matrix if null

    out <- list(
      Exposure = colnames_A,
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
    stop("The number of SNPs in betas.Gy, A matrix and the reference panel (Gl) are different.")
  }
}
