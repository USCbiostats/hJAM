#' Compute conditional Z matrix
#' @description The get_condZ function is to get the conditional Z matrix by using marginal Z matrix
#'
#' @param marginal_Z the marginal effects of SNPs on the exposures (Gx).
#' @param Gl the reference panel (Gl), such as 1000 Genome
#' @param N.Gx the sample size of each Gx. It can be a scalar or a vector. If there are multiple X's from different Gx, it should be a vector including the sample size of each Gx. If all alphas are from the same Gx, it could be a scalar.
#' @param ridgeTerm ridgeTerm = T when the matrix L is singular. Matrix L is obtained from the cholesky decomposition of G0'G0. Default as True.
#' @author Lai Jiang
#'
#' @examples
#' data(reference_data)
#' data(betas.Gy)
#' data(marginal_Z)
#' get_cond_Z(marginal_Z = marginal_Z, Gl = Gl, N.Gx = 1000, ridgeTerm = TRUE)

get_cond_Z =  function(marginal_Z, Gl, N.Gx, ridgeTerm = T){

  if(ncol(marginal_Z) == "NULL"){
    cat("! Please use get_cond_alpha instead of get_cond_Z.")
  }else if(N.Gx != ncol(marginal_Z) & length(N.Gx) != 1){
    cat("! ERROR: The length of the sample size of each Gx is different from the number of X in marginal Z matrix")
  }else if(nrow(marginal_Z) != ncol(Gl)){
    cat("! ERROR: The number of SNPs in marignal Z matrix and the reference panel (Gl) are different.")
  }else{

    # Check the length of N.Gx
    if(N.Gx != ncol(marginal_Z) & length(N.Gx) == 1){
      N.Gx = rep(N.Gx, ncol(marginal_Z))
    }

    # Check the dimension of Gl and Z
    dim_Gl = ncol(Gl)
    dim_Z = nrow(marginal_Z)

    if(dim_Gl == dim_Z){

      # Obtain the JAM variables: zL and L
      num_X = ncol(marginal_Z)
      conditional_Z = marginal_Z
      for(i_Z in 1:num_X){

        # Obtain marignal alpha and sample size for X_i
        alphas = marginal_Z[, i_Z]
        i_N = N.Gx[i_Z]

        # Compute the conditional alpha
        conditional_Z[, i_Z] = get_cond_alpha(alphas, Gl, N.Gx = i_N, ridgeTerm)
      }
      return(conditional_Z)
    }}
}
