#' Scatter plot for SNPs vs. one intermediate in the analysis
#' @description To generate the scatter plot of the SNPs vs. one intermediate that the user use in the analysis
#'
#' @param alphas The effects of SNPs on the intermediate (i.e. exposure/risk factor) (Gx).
#' @param betas.Gy The betas in the paper: the marginal effects of SNPs on the phenotype (Gy)
#' @param X.label The label of the intermediate (i.e. exposure/risk factor). Default is NULL.
#' @author Lai Jiang
#'
#' @return A set of scatter plots with x-axis being the conditional \eqn{\alpha}{alpha} estimates for each
#' intermediate and y-axis being the \eqn{\beta}{beta} estimates.
#' @export
#' @import ggplot2
#' @examples
#' data(MI)
#' t = SNPs_scatter_plot(alphas = MI.Amatrix[, 1], betas.Gy = MI.betas.gwas, X.label = "BMI")
#' t

SNPs_scatter_plot = function(alphas, betas.Gy, X.label = NULL){
  if(length(alphas) == length(betas.Gy)){
    dta.plot = data.frame(cbind(alphas, betas.Gy))
    p = ggplot(dta.plot, aes(x = alphas, y = betas.Gy)) + geom_point() +
      geom_hline(yintercept = 0, color = 'red') +
      xlab(expression(hat(alpha))) + ylab(expression(hat(beta))) + theme_bw()
    if(!is.null(X.label)){
      p = p + ggtitle(X.label)
    }
   return(p)
  }else{
    stop("The number of SNPs in alphas is not the same as the number of SNPs in betas.Gy vector.")
  }
}
