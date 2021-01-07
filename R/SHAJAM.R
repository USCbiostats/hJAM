#' SHA-JAM
#' Fit SHA-JAM
#' @description
#'
#' @param betas.Gy The betas in the paper: the marginal effects of SNPs on the phenotype (Gy)
#' @author Lai Jiang
#'
#' @return An object of the SHAJAM
#'
#' \describe{
#'    \item{Exposure}{The intermediates, such as the modifiable risk factors in Mendelian Randomization and gene expression in transcriptome analysis.}
#'    \item{numSNP}{The number of SNPs that the user use in the instrument set.}
#'    \item{Estimate}{The conditional estimates of the associations between intermediates and the outcome.}
#'    \item{StdErr}{The standard error of the conditional estimates of the associations between intermediates and the outcome.}
#'    \item{Lower.CI}{The lower bound of the 95\% confidence interval of the estimates.}
#'    \item{Upper.CI}{The upper bound of the 95\% confidence interval of the estimates.}
#'    \item{Pvalue}{The p value of the estimates with a type-I error equals 0.05.}
#' }
#'
#' @export
#' @import glmnet
#' @import susieR
#'

SHAJAM = function(betas.Gy, betas_se.Gy = NULL, N.Gy, MAF_summary = NULL,
                  Gl, A,
                  selection_alg = 'en',
                  L.susie = NULL, min_abs_corr = NULL, coverage=0.95,
                  estimate_residual_variance = TRUE,
                  max_iter = 500,
                  ridgeTerm = FALSE) {

  # Check the dimension of betas.Gy, Gl and A
  dim_betas = length(betas.Gy)
  dim_Gl = ncol(Gl)
  dim_A = ifelse(is.null(dim(A)), length(A), nrow(A))

  if(dim_betas == dim_Gl & dim_betas == dim_A){

    # Remove rows with all-zero
    zero.A.row = ifelse(apply(A, 1, sum)==0, TRUE, FALSE)
    A = A[!zero.A.row, ]
    betas.Gy = betas.Gy[!zero.A.row]
    betas_se.Gy = betas_se.Gy[!zero.A.row]
    Gl = Gl[, !zero.A.row]

    # Remove rows with zero in Genotype file
    if(sum(is.na(Gl))>0){
      Gl = Gl[complete.cases(Gl), ]
    }

    # Check the colnames of A matrix
    if(is.null(colnames(A))){
      stop("Please assign colnames to the A matrix input.\n")
    }

    if(is.null(MAF_summary)){
      p_D = apply(Gl, 2, mean)/2
    }else{
      p_D = MAF_summary[!zero.A.row]
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
    G0 = scale(Gl, center=T, scale=F)
    G0_t_G0 = t(G0)%*%G0

    ## Modify G0'G0 if the sample sizes of Gl and Gy are different
    Dj = 2*p_D*(1-p_D)*N.Gy
    D_sqrt = diag(sqrt(Dj))
    Dw_sqrt_inv = diag(1/sqrt(diag(G0_t_G0)))
    G0_t_G0.scaled = D_sqrt %*% Dw_sqrt_inv  %*% G0_t_G0 %*% Dw_sqrt_inv %*% D_sqrt

    Sj2 = betas_se.Gy^2
    yTy.est = median(Dj*Sj2*(N.Gy-1)+Dj*betas.Gy^2, na.rm = TRUE)

    if(selection_alg == 'susie'){
      if(is.null(L.susie) | is.null(min_abs_corr)){
        stop("Please specify the L value or min_abs_corr in Susie model.")
      }else if(is.null(betas_se.Gy)){
        stop("Specify the standard errors of the beta.gwas")
      }else{
        AtXtXA = t(A)%*%G0_t_G0.scaled%*%A
        XAty = t(A) %*% z

        susie_out = susieR::susie_suff_stat(XtX = AtXtXA, Xty = XAty,
                                            n = N.Gy, yty = yTy.est, L = L.susie,
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
          Selection_algorithm = selection_alg,
          num_Credible_sets = num.CS,
          all_variables = colnames(A),
          all_variable_pip = susie_out$pip,
          all_variable_coefficient = susieR::susie_get_posterior_mean(susie_out),
          numX = ncol(A),
          cs_purity = cs_purity)
      }
    }else{

      ## Add a ridge term in case G0'G0 is singular
      ridgeValue = ifelse(ridgeTerm, min(1, min(diag(G0_t_G0.scaled)*.001)), 0)
      G0_t_G0.ridge = G0_t_G0.scaled + ridgeValue*diag(length(betas.Gy))

      # Perfrom Cholesky decompostion and construct zL
      L = chol(G0_t_G0.ridge)
      zL = solve(t(L))%*%z

      # Perform linear regression
      X = L%*%A

      zero.X.column = ifelse(apply(X, 2, sum)==0, TRUE, FALSE)
      X = X[, !zero.X.column]
      if(selection_alg == 'en'){
        tune_alpha = 0.5
      }else if(selection_alg == 'lasso'){
        tune_alpha = 1
      }

      glm.out = cv.glmnet(X, zL, family="gaussian", alpha = tune_alpha, intercept = FALSE) # set intercept=0, use elastic net
      betas.XY = coef(glm.out)@x
      i.XY = coef(glm.out)@i
      i.length = length(i.XY)

      out <- list(
        numSNP = nrow(X),
        Selected_variable_length = i.length,
        Selected_variable_index = i.XY,
        Selected_variable_name = colnames(A)[i.XY],
        Coefficients = betas.XY,
        Selection_algorithm = selection_alg)
    }

    class(out) <- "SHAJAM"
    return(out)
  }else{
    stop("The number of SNPs in betas.Gy, A matrix and the reference panel (Gl) are different.")
  }
}
