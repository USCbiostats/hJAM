#' Elastic net hJAM
#' @description Function to implement regularized hJAM, including elastic net hJAM and lasso hJAM.
#'
#' @param betas.Gy The betas in the paper: the marginal effects of SNPs on the phenotype (Gy)
#' @param N.Gy The sample size of the GWAS where you obtain the betas.Gy and betas_se.Gy
#' @param eaf.Gy The effect allele frequency of the SNPs in betas.Gy
#' @param Geno The individual level data of the reference panel. Must have the same order of SNPs as in the betas.Gy.
#' @param A The conditional A matrix.
#' @param tune_glmnet The \eqn{\alpha} used in the glmnet R package to tune the shrinkage parameter. Default is 0.5.
#' @param ridgeTerm Add a small elelment to the diagnoal of X'X to make the matrix invertable.
#' @author Lai Jiang
#'
#' @return An object of the Regularized hJAM
#'
#' \describe{
#'    \item{numSNP}{The number of SNPs that the user use in the instrument set.}
#'    \item{Selected_variable_length}{The number of selected intermediates, regardless of the credible sets.}
#'    \item{Selected_variable_name}{The label/name for each selected intermediates.}
#'    \item{Coefficients}{The coefficients of selected intermediates. Otherwise will be zero.}
#' }
#'
#' @export
#' @importFrom glmnet cv.glmnet
#' @examples
#' data(ENhJAM.SimulationSet)
#' EN.hJAM(betas.Gy = Simulation.betas.gwas, N.Gy = 5000, eaf.Gy = Simulation.maf.gwas,
#' Geno = Simulation.Geno, A = Simulation.Amatrix, ridgeTerm = FALSE)

EN.hJAM = function(betas.Gy, N.Gy, eaf.Gy = NULL,
                   Geno, A, tune_glmnet = 0.5, ridgeTerm = FALSE){

  # Check the dimension of betas.Gy, Geno and A
  dim_betas = length(betas.Gy)
  dim_Geno = ncol(Geno)
  dim_A = ifelse(is.null(dim(A)), length(A), nrow(A))

  if(dim_betas == dim_Geno & dim_betas == dim_A){

    # Remove rows with all-zero
    zero.A.row = ifelse(apply(A, 1, sum)==0, TRUE, FALSE)
    A = A[!zero.A.row, ]
    betas.Gy = betas.Gy[!zero.A.row]
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

    glm.out = cv.glmnet(X, zL, family="gaussian", alpha = tune_glmnet, intercept = FALSE) # set intercept=0, use elastic net
    betas.XY = coef(glm.out)@x
    i.XY = coef(glm.out)@i
    i.length = length(i.XY)

    out <- list(
      numSNP = nrow(X),
      numX = ncol(X),
      Selected_variable_length = i.length,
      Selected_variable_index = i.XY,
      Selected_variable_name = colnames(A)[i.XY],
      Coefficients = betas.XY)

    class(out) <- "ENhJAM"
    return(out)
  }else{
    stop("The number of SNPs in betas.Gy, A matrix and the reference panel (Geno) are different.")
  }
}
