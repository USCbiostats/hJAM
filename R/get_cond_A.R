#' Compute conditional Z matrix
#' @description The get_cond_A function is to get the conditional A matrix by using marginal A matrix
#'
#' @param marginal_A the marginal effects of SNPs on the exposures (Gx).
#' @param Gl the reference panel (Gl), such as 1000 Genome
#' @param N.Gx the sample size of each Gx. It can be a scalar or a vector. If there are multiple X's from different Gx, it should be a vector including the sample size of each Gx. If all alphas are from the same Gx, it could be a scalar.
#' @param ridgeTerm ridgeTerm = TRUE when the matrix L is singular. Matrix L is obtained from the cholesky decomposition of G0'G0. Default as FALSE.
#' @author Lai Jiang
#'
#' @return A matrix with conditional estimates which are converted from marginal estimates using the JAM model.
#' @export
#' @examples
#' data(Gl)
#' data(betas.Gy)
#' data(marginal_A)
#' get_cond_A(marginal_A = marginal_A, Gl = Gl, N.Gx = c(339224, 659316), ridgeTerm = TRUE)

get_cond_A =  function(marginal_A, Gl, N.Gx, ridgeTerm = FALSE){

  if(ncol(marginal_A) == "NULL"){
    stop("Please use get_cond_alpha instead of get_cond_A.")
  }else if(length(N.Gx) != 1 && length(N.Gx) != ncol(marginal_A) ){
    stop("The length of the sample size of each Gx is different from the number of X in marginal A matrix")
  }else if(nrow(marginal_A) != ncol(Gl)){
    stop("The number of SNPs in marignal A matrix and the reference panel (Gl) are different.")
  }else{

    # Check the length of N.Gx
    if(length(N.Gx) != ncol(marginal_A) && length(N.Gx) == 1){
      N.Gx = rep(N.Gx, ncol(marginal_A))
    }

    # Check the dimension of Gl and A
    dim_Gl = ncol(Gl)
    dim_A = nrow(marginal_A)

    if(dim_Gl == dim_A){

      # Obtain the JAM variables: zL and L
      num_X = ncol(marginal_A)
      conditional_A = marginal_A
      for(i_A in 1:num_X){

        # Obtain marignal alpha and sample size for X_i
        alphas = marginal_A[, i_A]
        i_N = N.Gx[i_A]

        # Compute the conditional alpha
        conditional_A[, i_A] = get_cond_alpha(alphas, Gl, N.Gx = i_N, ridgeTerm)
      }
      return(conditional_A)
    }}
}

#' Compute conditional alphas
#' @description The get_cond_alpha function is to compute the conditional alpha vector for each X
#' If only one X in the model, please use get_cond_alpha instead of get_cond_A
#' A sub-step in the get_cond_A function
#'
#' @param alphas the marginal effects of SNPs on one exposure (Gx).
#' @param Gl the reference panel (Gl), such as 1000 Genome
#' @param N.Gx the sample size of the Gx. It can be a scalar.
#' @param ridgeTerm ridgeTerm = TRUE when the matrix L is singular. Matrix L is obtained from the cholesky decomposition of G0'G0. Default as FALSE
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
#' data(Gl)
#' data(betas.Gy)
#' data(marginal_A)
#' get_cond_alpha(alphas = marginal_A[, 1], Gl = Gl, N.Gx = 339224, ridgeTerm = TRUE)

get_cond_alpha = function(alphas, Gl, N.Gx, ridgeTerm = FALSE){

  ## Compute z vector
  p = apply(Gl, 2, mean)/2
  n0 = N.Gx*(1-p)^2
  n1 = N.Gx*2*p*(1-p)
  n2 = N.Gx*p^2

  y0 = -(n1*alphas+2*n2*alphas)/(n0+n1+n2)
  y1 = y0+alphas
  y2 = y0+2*alphas
  z = n1*y1 + 2*n2*y2

  ## Compute G0'G0
  G0 = scale(Gl, center=T, scale=F)
  G0_t_G0 = t(G0)%*%G0

  ## Modify G0'G0 if the sample sizes of Gl and Gx are different
  Dm = 2*p*(1-p)*N.Gx
  D_sqrt = diag(sqrt(Dm))
  Dw_sqrt_inv = diag(1/sqrt(diag(G0_t_G0)))
  G0_t_G0.scaled = D_sqrt %*% Dw_sqrt_inv  %*% t(G0) %*% G0 %*% Dw_sqrt_inv %*% D_sqrt

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
