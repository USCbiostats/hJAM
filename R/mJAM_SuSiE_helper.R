#' Get SuSiE posterior mean
#'
#' @param res A SuSiE fit object
#' @param prior_tol When the prior variance is estimated, compare the estimated value to prior_tol at the end of the computation, and exclude a single effect from PIP computation if the estimated prior variance is smaller than this tolerance value.
#'
#' @returns A vector of posterior mean effects

susie_get_posterior_mean_v2 <- function (res, prior_tol = 1e-09)
{
  if (is.numeric(res$V)) {
    include_idx = which(res$V > prior_tol)
  }
  else {include_idx = 1:nrow(res$alpha)}

  if (length(include_idx) > 1) {
    return(colSums((res$alpha * res$mu)[include_idx, ])/res$X_column_scale_factors)
  } else if (length(include_idx) == 1) {
    return((res$alpha * res$mu)[include_idx, ]/res$X_column_scale_factors)
  }else {
    return(numeric(ncol(res$mu)))
  }
}

#' Get SuSiE posterior sd
#'
#' @param res A SuSiE fit object
#' @param prior_tol When the prior variance is estimated, compare the estimated value to prior_tol at the end of the computation, and exclude a single effect from PIP computation if the estimated prior variance is smaller than this tolerance value.
#'
#' @returns A vector of posterior standard deviations

susie_get_posterior_sd_v2 <- function (res, prior_tol = 1e-09)
{
  if (is.numeric(res$V))
    include_idx = which(res$V > prior_tol)
  else include_idx = 1:nrow(res$alpha)
  if (length(include_idx) > 1) {
    return(sqrt(colSums((res$alpha * res$mu2 - (res$alpha *
                                                  res$mu)^2)[include_idx, ]))/(res$X_column_scale_factors))
  } else if (length(include_idx) == 1){
    return(sqrt((res$alpha * res$mu2 - (res$alpha * res$mu)^2)[include_idx, ])/(res$X_column_scale_factors))
  } else {return(numeric(ncol(res$mu)))}
}

#' Get and tidy SuSiE credible sets
#'
#' @param susie_fit A SuSiE fit object
#' @param SNP_names A character vector of the names of the SNPs. Note that the order of SNP names should be the same as `EAFs`, `X_ref`, `marginal.betas` and `marginal.se` in `mJAM_SuSiE()`.
#' @param coverage A number between 0 and 1 specifying the “coverage” of the estimated confidence sets.
#'
#' @author Jiayi Shen
#' @export
#'
#' @import dplyr
#' @returns A table summary of SuSiE credible sets with the following columns:
#'
#'#' \describe{
#'    \item{index}{The label for a distinct credible set.}
#'    \item{coverage}{The empirical coverage of this credible set.}
#'    \item{CS_size}{The number of SNPs in total in corresponding credible set.}
#'    \item{index_SNP_id}{The name of the index SNP (SNP with highest posterior probability) in corresponding credible set.}
#'    \item{CS_SNP_id}{The names of individual SNPs selected in this credible set.}
#' }

mJAM_SuSiE_get_cs <- function(susie_fit,susie_pip,
                              coverage = 0.95){

  SNP_names <- rownames(susie_pip)

  # cs_output <- susieR::susie_get_cs(susie_fit, coverage = coverage)
  cs_output <- susie_fit$sets
  if (is.null(cs_output$cs)){
    cs_summary <- data.frame(index = NA,
                             coverage = coverage,
                             CS_size = 0,
                             index_SNP = "None Selected",
                             CS_SNP = "None Selected"
    )
  }else{
    cs_summary <- data.frame(index = character(),coverage = double(),
                             CS_size = integer(),index_SNP_id = integer(),
                             CS_SNP_id = integer())
    for(r in 1:length(cs_output$cs)){
      ## identify index SNP
      temp_index_id <- rownames(susie_pip[cs_output$cs[[r]],])[which.max(susie_pip[cs_output$cs[[r]],]$pip)]
      ## pull CS SNPs
      temp_cs_summary <- data.frame(index = names(cs_output$cs)[r],
                                    coverage = cs_output$coverage[r],
                                    CS_size = length(cs_output$cs[[r]]),
                                    index_SNP = temp_index_id,
                                    CS_SNP_id = cs_output$cs[[r]]
      )
      ## concat
      cs_summary <- rbind(cs_summary,temp_cs_summary)
    }
    cs_summary <- cs_summary %>%
      left_join(data.frame(CS_SNP_id = 1:length(SNP_names), CS_SNP = SNP_names), by = "CS_SNP_id") %>%
      dplyr::select(-c(CS_SNP_id))

  }

  return(cs_summary)

}


