#' Heatmap for all the SNPs used in the analysis
#' @description To generate the heatmap of all the SNPs that the user use in the analysis
#'
#' @param Gl The reference panel (Gl) of the SNPs that the user use in the analysis, such as 1000 Genome
#' @author Lai Jiang
#'
#' @export
#' @import dplyr reshape2
#' @importFrom stats cor
#' @examples
#' data(Gl)
#' t = SNPs_heatmap(Gl = Gl)
#' t

SNPs_heatmap = function(Gl){
  rho = cor(Gl)
  melted_cormat = melt(rho)
  heatmap_p = ggplot(data = melted_cormat, aes(melted_cormat[, 1], melted_cormat[, 2], fill = melted_cormat[, 3]))+
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal()+
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
    coord_fixed() + xlab(" ") + ylab(" ")
  return(heatmap_p)
}
