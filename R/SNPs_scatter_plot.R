#' Scatter plot for all the SNPs used in the analysis
#' @description To generate the scatter plot of all the SNPs that the user use in the analysis
#'
#' @param A The effects of SNPs on the exposures (Gx).
#' @param betas.Gy The betas in the paper: the marginal effects of SNPs on the phenotype (Gy)
#' @param num_X The number of intermediates in the research question.
#' @author Lai Jiang
#'
#' @return A set of scatter plots with x-axis being the conditional \eqn{\alpha}{alpha} estimates for each
#' intermediate and y-axis being the \eqn{\beta}{beta} estimates.
#' @export
#' @import ggplot2 ggpubr
#' @examples
#' data(conditional_A)
#' data(betas.Gy)
#' t = SNPs_scatter_plot(A = conditional_A, betas.Gy = betas.Gy, num_X = 2)
#' t

SNPs_scatter_plot = function(A, betas.Gy, num_X){
  plotList = list()
  label.plot = colnames(A)
  if(nrow(A) == length(betas.Gy)){
    dta.plot = as.data.frame(cbind(A, betas.Gy))
    for(i in 1:num_X){
      plotList[[i]] = ggplot(dta.plot, aes(x = dta.plot[, i], y = betas.Gy)) + geom_point() +
                     geom_hline(yintercept = 0, color = 'red') +
                     xlab(expression(hat(alpha))) + ylab(expression(hat(beta))) + theme_bw()
   }
   p = ggarrange(plotlist=plotList, ncol=2, nrow=ceiling(num_X/2), common.legend = T, legend="bottom",
                  labels = label.plot, font.label = list(size = 12))
   return(p)
  }else{
    stop("The number of SNPs in A matrix is not the same as the number of SNPs in betas.Gy vector.")
  }
}
