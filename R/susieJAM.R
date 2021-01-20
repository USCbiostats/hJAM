#' Compute conditional A using SuSiE JAM
#'
#' @description The \code{susieJAM_A} function is to get the conditional A matrix by using marginal A matrix
#'
#' @param marginalA the marginal effects of SNPs on the exposures (Gx).
#' @param marginalA_se the standard error of the marginal effects of SNPs on the exposures (Gx).
#' @param N.Gx the sample size of each Gx. It can be a scalar or a vector. If there are multiple X's from different Gx, it should be a vector including the sample size of each Gx. If all alphas are from the same Gx, it could be a scalar.
#' @param raf.Gy the effect allele frequency of the SNPs in the Gx data.
#' @param Geno the reference panel (Geno), such as 1000 Genome。
#' @param inclusion.indicator The matrix of inclusion indicator of SNPs for each intermediate. Included as 1; otherwise 0.
#' @param L.cs A susie input parameter. Number of components (nonzero elements) in the SuSiE regression model. If L.cs is larger than the number of covariate (p), L.cs is set to p.
#' @param min_abs_corr A susie input parameter. Minimum of absolute value of correlation allowed in a credible set. The default, 0.5, corresponds to squared correlation of 0.25, which is a commonly used threshold for genotype data in genetics studies.
#' @param max_iter Maximum number of iterations in SuSiE fitting.
#' @param coverage Default as 0.95.The coveralge level of the credible set.
#' @param estimate_residual_variance Default as TRUE. Estimate the residual variance in each iteration of SuSiE fitting.
#'
#' @author Lai Jiang
#'
#' @return A matrix with conditional estimates which are converted from marginal estimates using the susie JAM model.
#' @export
#' @importFrom stats sd
#' @examples
#' data(GTEx.PrCa)
#' susieJAM_A(marginalA = GTEx.PrCa.marginal.A[, 1:10],
#' marginalA_se = GTEx.PrCa.marginal.A.se[, 1:10], raf.Gy = GTEx.PrCa.maf.gwas,
#' Geno = GTEx.PrCa.Geno, inclusion.indicator = GTEx.PrCa.inclusion.indicator,
#' N.Gx = 620, L.cs = 10, min_abs_corr = 0.5)

susieJAM_A = function(marginalA, marginalA_se, N.Gx, raf.Gy = NULL, Geno,
                      inclusion.indicator,
                      L.cs, min_abs_corr, max_iter, coverage,
                      estimate_residual_variance=TRUE){

  if(ncol(marginalA) == "NULL"){
    stop("Please use susieJAM_alphas instead of susieJAM_A")
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

      num_X = ncol(marginalA)
      cond_A = marginalA

      for(i_A in 1:ncol(marginalA)){

        i_row = which(inclusion.indicator[, i_A] == 1)
        alphas = marginalA[i_row, i_A]
        alphas_se = marginalA_se[i_row, i_A]
        i_raf = raf.Gy[i_row]
        i_Geno = Geno[, i_row]

        if(length(alphas) > 1){
          i_alpha = susieJAM_alphas(marginalA = alphas, marginalA_se = alphas_se, raf.Gy = i_raf,
                                    Geno = i_Geno, N.Gx = 620, L.cs = 10, min_abs_corr = 0.5)
          cond_A[i_row, i_A] = i_alpha$cond_alphas
        }else{
          cond_A[i_row, i_A] = alphas
        }
      }

      return(cond_A)
    }
  }
}

#' Compute conditional alphas using SuSiE JAM
#' @description The \code{susieJAM_alphas} function is to perform the variable selection and compute the selected conditional alpha vector for one intermediate.
#' If only one intermediate in the model, please use susieJAM_alphas instead of susieJAM_A
#'
#' @param marginalA the marginal effects of SNPs on one exposure (Gx).
#' @param marginalA_se the standard error of the marginal effects of SNPs on one outcome (Gx).
#' @param Geno the reference panel (Geno), such as 1000 Genome. The reference data has to be centered.
#' @param N.Gx the sample size of the Gx. It can be a scalar.
#' @param raf.Gy The vector of the minor allele frequency or effect allele frequency in the GWAS.
#' @param L.cs A susie input parameter. Number of components (nonzero elements) in the SuSiE regression model. If L.cs is larger than the number of covariate (p), L.cs is set to p.
#' @param min_abs_corr A susie input parameter. Minimum of absolute value of correlation allowed in a credible set. The default, 0.5, corresponds to squared correlation of 0.25, which is a commonly used threshold for genotype data in genetics studies.
#' @param max_iter Maximum number of iterations in SuSiE fitting.
#' @param coverage Default as 0.95.The coveralge level of the credible set.
#' @param estimate_residual_variance Default as TRUE. Estimate the residual variance in each iteration of SuSiE fitting.
#'
#' @author Lai Jiang
#'
#' @export
#' @examples
#' data(GTEx.PrCa)
#' include.SNPs = which(GTEx.PrCa.inclusion.indicator[,1]==1)
#' susieJAM_alphas(marginalA = GTEx.PrCa.marginal.A[include.SNPs, 1],
#' marginalA_se = GTEx.PrCa.marginal.A.se[include.SNPs, 1], raf.Gy = GTEx.PrCa.maf.gwas[include.SNPs],
#' Geno = GTEx.PrCa.Geno[, include.SNPs], N.Gx = 620, L.cs = 10, min_abs_corr = 0.5)

susieJAM_alphas = function(marginalA, marginalA_se, N.Gx, raf.Gy = NULL, Geno,
                           L.cs = 10, min_abs_corr = 0.6, max_iter = 100,
                           coverage = 0.95,
                           estimate_residual_variance = FALSE){

  # Define alphas and alphas_se
  alphas = marginalA
  alphas_se = marginalA_se

  # Check the reference panel component
  if(is.null(Geno)){
    stop("ERROR: Need Geno or the correlation coefficient matrix of centered Geno")
  }

  ## Remove NA in Geno
  if(sum(is.na(Geno))>0){Geno = Geno[complete.cases(Geno), ]}

  # Remove rows without variation or missing data
  i_sd0 = ifelse(apply(Geno, 2, sd) != 0, TRUE, FALSE)
  Geno = Geno[, i_sd0]
  alphas = alphas[i_sd0]
  alphas_se = alphas_se[i_sd0]
  i_raf = raf.Gy[i_sd0]

  # Get the JAM components
  p_D = raf.Gy
  n0 = N.Gx*(1-p_D)^2
  n1 = N.Gx*2*p_D*(1-p_D)
  n2 = N.Gx*p_D^2

  y0 = -(n1*alphas+2*n2*alphas)/(n0+n1+n2)
  y1 = y0+alphas
  y2 = y0+2*alphas
  z = n1*y1 + 2*n2*y2

  ## Modify G0'G0 if the sample sizes of Geno and Gx are different
  Geno = Geno[complete.cases(Geno), ]
  Dj = 2*p_D*(1-p_D)*N.Gx
  D_sqrt = diag(sqrt(Dj))
  rho_Geno = WGCNA::cor(Geno)
  G0_t_G0 = D_sqrt %*% rho_Geno %*% D_sqrt

  Sj2 = alphas_se^2
  yTy.est = median(Dj*Sj2*(N.Gx-1)+Dj*alphas^2)

  susie_out = susieR::susie_suff_stat(XtX = G0_t_G0, Xty = z,
                                      n = N.Gx, yty = yTy.est, L = L.cs,
                                      min_abs_corr = min_abs_corr, max_iter = max_iter,
                                      estimate_residual_variance = estimate_residual_variance,
                                      coverage = coverage)
  cond_alphas = susie_get_posterior_mean(susie_out)

  return(list(susie.out = susie_out,
              cond_alphas = cond_alphas,
              pip = susie_out$pip))
}
