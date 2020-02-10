#' Keep the output as three digits
#' @description Keep the output as three digits
#' @author Lai Jiang
#'
#' @param x input
#' @param ... other options you want to put in
#' @export
output.format = function(x, ...){
  return(round(x, digits = 3))
}

#' Print out for hJAM_lnreg
#' @description Print out for hJAM_lnreg
#' @author Lai Jiang
#'
#' @param x input
#' @param ... other options you want to put in
#' @export
print.hJAM_lnreg = function(x, ...) {
  cat("------------------------------------------------------", "\n")
  cat("                   hJAM output                        ", "\n")
  cat("------------------------------------------------------", "\n")

  cat(paste0("Number of SNPs used in model: ", x$numSNP), "\n\n")

  output =  as.data.frame(unclass(x))[c(3:7)]
  output$'95% CI' = paste0("(", output.format(output$Lower.CI), ", ", output.format(output$Upper.CI), ")")

  outprint = output[, c('Estimate', 'StdErr', '95% CI', 'Pvalue')]
  outprint[, c(1:2)] = sapply(outprint[, c(1,2)], function(x) output.format(x))
  rownames(outprint) <- x$Exposure

  print(outprint)
  cat("------------------------------------------------------", "\n\n")
}

#' Print out for hJAM_egger
#' @description Print out for hJAM_egger
#' @author Lai Jiang
#'
#' @param x input
#' @param ... other options you want to put in
#' @export
print.hJAM_egger = function(x, ...) {
  cat("------------------------------------------------------", "\n")
  cat("                   hJAM egger output                  ", "\n")
  cat("------------------------------------------------------", "\n")

  cat(paste0("Number of SNPs used in model: ", x$numSNP), "\n\n")

  output =  as.data.frame(unclass(x))[c(3:12)]
  output$'95% CI' = paste0("(", output.format(output$Lower.CI), ", ", output.format(output$Upper.CI), ")")
  int.95ci = paste0("(", output.format(x$Lower.CI.Int), ", ", output.format(x$Upper.CI.Int), ")")

  exp_output = output[, c('Estimate', 'StdErr', '95% CI', 'Pvalue')]
  exp_output[, c(1:2)] = sapply(exp_output[, c(1,2)], function(x) output.format(x))
  rownames(exp_output) <- x$Exposure
  print(exp_output)

  cat("\nIntercept\n")
  int_output = cbind(as.numeric(output.format(x$Est.Int)), as.numeric(output.format(x$StdErr.Int)),
                     int.95ci, as.numeric(output.format(x$Pvalue.Int)))
  colnames(int_output) = c("Est.Int", "StdErr.Int", "95% CI.Int", "Pvalue.Int")
  print(int_output)
  cat("------------------------------------------------------", "\n\n")
}
