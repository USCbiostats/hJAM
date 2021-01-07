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

#' Print out for hJAM
#' @description Print out for hJAM_lnreg
#' @author Lai Jiang
#'
#' @param x object output by hJAM
#' @param ... other options you want to put in
#' @export
print.hJAM = function(x, ...) {
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
#' @param x obejct output from hJAM_egger
#' @param ... other options you want to put in
#' @export
print.hJAM_egger = function(x, ...) {
  cat("------------------------------------------------------", "\n")
  cat("                   hJAM egger output                  ", "\n")
  cat("------------------------------------------------------", "\n")

  cat(paste0("Number of SNPs used in model: ", x$numSNP), "\n\n")

  output =  as.data.frame(unclass(x)[c(3:7)])
  output$'95% CI' = paste0("(", output.format(output$Lower.CI), ", ", output.format(output$Upper.CI), ")")
  int.95ci = paste0("(", output.format(x$Lower.CI.Int), ", ", output.format(x$Upper.CI.Int), ")")

  exp_output = output[, c('Estimate', 'StdErr', '95% CI', 'Pvalue')]
  exp_output[, c(1:2, 4)] = sapply(exp_output[, c(1,2, 4)], function(x) output.format(x))

  int_output = cbind(as.numeric(output.format(x$Est.Int)), as.numeric(output.format(x$StdErr.Int)),
                     int.95ci, as.numeric(output.format(x$Pvalue.Int)))
  colnames(int_output) = colnames(exp_output)
  exp_output = rbind(exp_output, int_output)
  rownames(exp_output) <- c(x$Exposure, 'Intercept')
  print(exp_output)

  cat("------------------------------------------------------", "\n\n")
}

#' Print out for SHA-JAM
#' @description Print out for SHA-JAM
#' @author Lai Jiang
#'
#' @param x obejct output from SHAJAM
#' @param ... other options you want to put in
#' @export
print.SHAJAM = function(x, ...) {
  cat("------------------------------------------------------", "\n")
  cat("                   SHAJAM output                      ", "\n")
  cat("------------------------------------------------------", "\n")

  cat(paste0("Number of SNPs used in model: ", x$numSNP), "\n")
  cat(paste0("Number of intermediates used in model: ", x$numX), "\n")

  if(x$Selection_algorithm == 'susie'){
    cat(paste0("Number of the credible sets of intermediates selected by ", x$Selection_algorithm, ": ", x$num_Credible_sets), "\n\n")

    if(x$Selected_variable_length == 0){
      cat(paste0("No variable has been selected by ", x$Selection_algorithm, "\n"))
    }else{
      output =  as.data.frame(unclass(x)[4:7])
      output[, 3:4] = sapply(output[, 3:4], function(x) output.format(x))
      colnames(output) = c('Credible Sets', 'Variable', 'Coefficients', 'PIP')
      rownames(output) = NULL
      print(output)
    }
  }else{
    cat(paste0("Number of intermediates selected by ", x$Selection_algorithm, ": ", x$Selected_variable_length), "\n\n")

    if(x$Selected_variable_length == 0){
      cat(paste0("No variable has been selected by ", x$Selection_algorithm, "\n"))
    }else{
      output =  as.data.frame(unclass(x))[4:5]
      output[, 2] = sapply(output[, 2], function(x) output.format(x))
      colnames(output) = c('Variable', 'Coefficients')
      print(output)
    }
  }

  cat("------------------------------------------------------", "\n\n")
}
