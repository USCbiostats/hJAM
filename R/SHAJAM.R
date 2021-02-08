#' SHA-JAM
#' Fit SHA-JAM
#' @description Function to implement SHA-JAM
#'
#' @param betas.Gy The betas in the paper: the marginal effects of SNPs on the phenotype (Gy)
#' @param betas_se.Gy The standard errors of the betas
#' @param N.Gy The sample size of the GWAS where you obtain the betas.Gy and betas_se.Gy
#' @param eaf.Gy The reference allele frequency of the SNPs in betas.Gy
#' @param Geno The individual level data of the reference panel. Must have the same order of SNPs as in the betas.Gy.
#' @param A The conditional A matrix.
#' @param L.cs The largest number of credible set allowed in SHA-JAM. Required by SHA-JAM.
#' @param min_abs_corr The requested minimum absolute correlation coefficient between intermediates within one credible set. Required by SHA-JAM.
#' @param coverage The coverage of credible set. Default is 0.95. Required by SHA-JAM.
#' @param estimate_residual_variance If estimate the residual variance in the fitting procedure of SHA-JAM. Default as TRUE. Required by SHA-JAM.
#' @param max_iter The number of maximum iterations in fitting SHA-JAM. Required by SHA-JAM.
#' @author Lai Jiang
#'
#' @return An object of the SHAJAM
#'
#' \describe{
#'    \item{numSNP}{The number of SNPs used in the analysis.}
#'    \item{numX}{The number of intermediates in the analysis.}
#'    \item{Selected_variable_length}{The number of selected intermediates, regardless of the credible sets.}
#'    \item{Selected_variable_name}{The label/name for each selected intermediates.}
#'    \item{Coefficients}{The coefficients of selected intermediates.}
#'    \item{Selected_variable_pip}{The posterior inclusion probability of each selected intermediate.}
#'    \item{num_Credible_sets}{Number of credible sets.}
#'    \item{all_variables}{The label/name for all candidate intermediates.}
#'    \item{all_variable_pip}{The posterior inclusion probability of all candidate intermediates.}
#'    \item{all_variable_coefficient}{The coefficients of all candidate intermediates.}
#'    \item{cs_purity}{The purity of the credibel set selected.}
#' }
#'
#' @export
#' @importFrom stats coef complete.cases median
#' @import susieR

SHAJAM = function(betas.Gy, betas_se.Gy = NULL, N.Gy, eaf.Gy = NULL,
                  Geno, A, L.cs = NULL, min_abs_corr = NULL, coverage=0.95,
                  estimate_residual_variance = TRUE,
                  max_iter = 500) {

  # Check the dimension of betas.Gy, Geno and A
  dim_betas = length(betas.Gy)
  dim_Geno = ncol(Geno)
  dim_A = ifelse(is.null(dim(A)), length(A), nrow(A))

  if(dim_betas == dim_Geno & dim_betas == dim_A){

    # Remove rows with all-zero
    zero.A.row = ifelse(apply(A, 1, sum)==0, TRUE, FALSE)
    A = A[!zero.A.row, ]
    betas.Gy = betas.Gy[!zero.A.row]
    betas_se.Gy = betas_se.Gy[!zero.A.row]
    Geno = Geno[, !zero.A.row]

    # Remove rows with zero in Genotype file
    if(sum(is.na(Geno))>0){
      Geno = Geno[complete.cases(Geno), ]
    }

    # Check the colnames of A matrix
    if(is.null(colnames(A))){
      stop("Please assign colnames to the A matrix input.\n")
    }

    if(is.null(eaf.Gy)){
      p_D = apply(Geno, 2, mean)/2
    }else{
      p_D = eaf.Gy[!zero.A.row]
    }

    # Obtain the JAM variables: zL and L
    n0 = N.Gy*(1-p_D)^2
    n1 = N.Gy*2*p_D*(1-p_D)
    n2 = N.Gy*p_D^2

    y0 = -(n1*betas.Gy+2*n2*betas.Gy)/(n0+n1+n2)
    y1 = y0+betas.Gy
    y2 = y0+2*betas.Gy
    z = n1*y1 + 2*n2*y2

    ## Compute G0'G0
    G0 = scale(Geno, center=T, scale=F)
    G0_t_G0 = t(G0)%*%G0

    ## Modify G0'G0 if the sample sizes of Geno and Gy are different
    Dj = 2*p_D*(1-p_D)*N.Gy
    D_sqrt = diag(sqrt(Dj))
    Dw_sqrt_inv = diag(1/sqrt(diag(G0_t_G0)))
    G0_t_G0.scaled = D_sqrt %*% Dw_sqrt_inv  %*% G0_t_G0 %*% Dw_sqrt_inv %*% D_sqrt

    Sj2 = betas_se.Gy^2
    yTy.est = median(Dj*Sj2*(N.Gy-1)+Dj*betas.Gy^2, na.rm = TRUE)

    if(is.null(L.cs) | is.null(min_abs_corr)){
      stop("Please specify the L value or min_abs_corr in Susie model.")
    }else if(is.null(betas_se.Gy)){
      stop("Specify the standard errors of the beta.gwas")
    }else{
      AtXtXA = t(A)%*%G0_t_G0.scaled%*%A
      XAty = t(A) %*% z

      susie_out = susieR::susie_suff_stat(XtX = AtXtXA, Xty = XAty,
                                          n = N.Gy, yty = yTy.est, L = L.cs,
                                          min_abs_corr = min_abs_corr, max_iter = max_iter, coverage=coverage,
                                          estimate_residual_variance = estimate_residual_variance)
      num.CS = length(susie_out$sets$cs)
      cs_purity = susie_out$sets$purity

      if(num.CS > 0){
        i.CS.name = names(eval(parse(text = "susie_out$sets$cs[1]")))
        Selected_credible_sets = rep(i.CS.name, times = eval(parse(text = paste0("length(susie_out$sets$cs$", i.CS.name, ")"))))
        i.XY = unlist(eval(parse(text = "susie_out$sets$cs[1]")))
        if(num.CS > 1){
          for(j.cs in 2:num.CS){
            i.CS.name = names(eval(parse(text = paste0("susie_out$sets$cs[", j.cs, "]"))))
            Selected_credible_sets = c(Selected_credible_sets,
                                       rep(i.CS.name, times = eval(parse(text = paste0("length(susie_out$sets$cs$", i.CS.name, ")")))))
            i.XY = c(i.XY, unlist(eval(parse(text = paste0("susie_out$sets$cs[", j.cs, "]")))))
          }
        }
      }else{
        Selected_credible_sets = i.XY = NULL
      }

      betas.XY = susieR::susie_get_posterior_mean(susie_out)[i.XY]
      i.length = length(unique(unlist(susie_out$sets$cs)))
      i.pip = susie_out$pip[i.XY]

      out <- list(
        numSNP = length(betas.Gy),
        Selected_variable_length = i.length,
        Selected_variable_index = i.XY,
        Selected_credible_sets = Selected_credible_sets,
        Selected_variable_name = colnames(A)[i.XY],
        Coefficients = betas.XY,
        Selected_variable_pip = i.pip,
        num_Credible_sets = num.CS,
        all_variables = colnames(A),
        all_variable_pip = susie_out$pip,
        all_variable_coefficient = susieR::susie_get_posterior_mean(susie_out),
        numX = ncol(A),
        cs_purity = cs_purity)
    }

    class(out) <- "SHAJAM"
    return(out)
  }else{
    stop("The number of SNPs in betas.Gy, A matrix and the reference panel (Geno) are different.")
  }
}
