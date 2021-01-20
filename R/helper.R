#' Heatmap for all the SNPs used in the analysis
#' @description To generate the heatmap of all the SNPs that the user use in the analysis
#'
#' @param Geno The reference panel (Geno) of the SNPs that the user use in the analysis, such as 1000 Genome
#' @param show.variables Select to show the variables name or not. Default set to be FALSE.
#' @author Lai Jiang
#'
#' @export
#' @importFrom reshape2 melt
#' @importFrom WGCNA cor
#' @examples
#' data(MI.Rdata)
#' t = SNPs_heatmap(Geno = MI.Geno[, 1: 10])
#' t
#' t = SNPs_heatmap(Geno = MI.Geno[, 1: 10], show.variable = TRUE)
#' t

SNPs_heatmap = function(Geno, show.variables = FALSE){
  rho = WGCNA::cor(Geno)
  melted_cormat = reshape2::melt(rho)
  melted_cormat$r.text = paste0("r = ", round(melted_cormat$value, 3))

  heatmap_p = ggplot(data = melted_cormat, aes(y=melted_cormat[, 1], x=melted_cormat[, 2], fill = melted_cormat[, 3]))+
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90)) +
    coord_fixed() + xlab(" ") + ylab(" ")

  if(!show.variables){
    heatmap_p = heatmap_p + theme(axis.text = element_blank(), axis.line=element_blank(),
                                  axis.ticks = element_blank())
  }
  heatmap_p
  return(heatmap_p)
}

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
