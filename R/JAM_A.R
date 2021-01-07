#' Compute conditional A matrix
#' @description The \code{JAM_A} function is to get the conditional A matrix by using marginal A matrix
#'
#' @param marginalA the marginal effects of SNPs on the exposures (Gx).
#' @param Geno the reference panel (Geno), such as 1000 Genome
#' @param N.Gx the sample size of each Gx. It can be a scalar or a vector. If there are multiple X's from different Gx, it should be a vector including the sample size of each Gx. If all alphas are from the same Gx, it could be a scalar.
#' @param eaf_Gx the effect allele frequency of the SNPs in the Gx data.
#' @param ridgeTerm ridgeTerm = TRUE when the matrix L is singular. Matrix L is obtained from the cholesky decomposition of G0'G0. Default as FALSE.
#' @author Lai Jiang
#'
#' @return A matrix with conditional estimates which are converted from marginal estimates using the JAM model.
#' @export
#' @examples
#' data(MI)
#' JAM_A(marginalA = MI.marginal.Amatrix, Geno = MI.Geno, N.Gx = c(339224, 659316), ridgeTerm = TRUE)
#' JAM_A(marginalA = MI.marginal.Amatrix, Geno = MI.Geno, N.Gx = c(339224, 659316), eaf_Gx = MI.SNPs_info$ref_frq, ridgeTerm = TRUE)

JAM_A =  function(marginalA, Geno, N.Gx, eaf_Gx = NULL, ridgeTerm = FALSE){

  if(ncol(marginalA) == "NULL"){
    stop("Please use JAM_alphas instead of get_cond_A.")
  }else if(length(N.Gx) != 1 && length(N.Gx) != ncol(marginalA) ){
    stop("The length of the sample size of each Gx is different from the number of X in marginal A matrix")
  }else if(nrow(marginalA) != ncol(Geno)){
    stop("The number of SNPs in marignal A matrix and the reference panel (Geno) are different.")
  }else{

    # Check the length of N.Gx
    if(length(N.Gx) != ncol(marginalA) && length(N.Gx) == 1){
      N.Gx = rep(N.Gx, ncol(marginalA))
    }

    # Check the dimension of Geno and A
    dim_Geno = ncol(Geno)
    dim_A = nrow(marginalA)

    if(dim_Geno == dim_A){

      # Obtain the JAM variables: zL and L
      num_X = ncol(marginalA)
      conditional_A = marginalA
      for(i_A in 1:num_X){

        # Obtain marignal alpha and sample size for X_i
        alphas = marginalA[, i_A]
        i_N = N.Gx[i_A]

        # Compute the conditional alpha
        conditional_A[, i_A] = JAM_alphas(alphas, Geno, N.Gx = i_N, eaf_Gx = eaf_Gx, ridgeTerm)
      }
      return(conditional_A)
    }
    }
}

#' Compute conditional alphas
#' @description The \code{JAM_alphas} function is to compute the conditional alpha vector for each X
#' If only one X in the model, please use JAM_alphas instead of JAM_A
#' A sub-step in the JAM_A function
#'
#' @param alphas the marginal effects of SNPs on one exposure (Gx).
#' @param Geno the reference panel (Geno), such as 1000 Genome
#' @param N.Gx the sample size of the Gx. It can be a scalar.
#' @param eaf_Gx the effect allele frequency of the SNPs in the Gx data.
#' @param ridgeTerm ridgeTerm = TRUE when the matrix L is singular. Matrix L is obtained from the cholesky decomposition of G0'G0. Default as FALSE
#' @param eaf_Gx TBD
#' @author Lai Jiang
#' @return A vector with conditional estimates which are converted from marginal estimates using the JAM model.
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
#' data(MI)
#' JAM_alphas(alphas = MI.marginal.Amatrix[, 1], Geno = MI.Geno, N.Gx = 339224, ridgeTerm = TRUE)
#' JAM_alphas(alphas = MI.marginal.Amatrix[, 1], Geno = MI.Geno, N.Gx = 339224, eaf_Gx = MI.SNPs_info$ref_frq, ridgeTerm = TRUE)

JAM_alphas = function(alphas, Geno, N.Gx, eaf_Gx = NULL, ridgeTerm = FALSE){

  ## Compute z vector
  if(!is.null(eaf_Gx)){
    p = eaf_Gx
  }else{
    p = apply(Geno, 2, mean)/2
  }

  n0 = N.Gx*(1-p)^2
  n1 = N.Gx*2*p*(1-p)
  n2 = N.Gx*p^2

  y0 = -(n1*alphas+2*n2*alphas)/(n0+n1+n2)
  y1 = y0+alphas
  y2 = y0+2*alphas
  z = n1*y1 + 2*n2*y2

  ## Compute G0'G0
  G0 = scale(Geno, center=T, scale=F)
  G0_t_G0 = t(G0)%*%G0

  ## Modify G0'G0 if the sample sizes of Geno and Gx are different
  Dm = 2*p*(1-p)*N.Gx
  D_sqrt = diag(sqrt(Dm))
  Dw_sqrt_inv = diag(1/sqrt(diag(G0_t_G0)))
  G0_t_G0.scaled = D_sqrt %*% Dw_sqrt_inv  %*% G0_t_G0 %*% Dw_sqrt_inv %*% D_sqrt

  ## Add a ridge term in case G0'G0 is singular
  ridgeValue = ifelse(ridgeTerm, min(1,min(diag(G0_t_G0.scaled)*.001)), 0)
  G0_t_G0.ridge = G0_t_G0.scaled + ridgeValue*diag(length(alphas))

  ## Get L matrix and zL vector
  L = chol(G0_t_G0.ridge)
  zL = solve(t(L))%*%z

  # Get the conditional alpha vectors
  cond_alphas = summary(stats::lm(zL ~ 0 + L))$coef[ ,1]

  return(cond_alphas)
}
