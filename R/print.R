#' @export
print.hJAM_lnreg <- function(x, ...) {
  cat("--------------------------------------------", "\n")
  cat("                hJAM output                 ", "\n")
  cat("--------------------------------------------", "\n")

  cat(paste0("Number of SNPs used in model: ", x$numSNP), "\n\n")

  output =  as.data.frame(unclass(x))[c(3, 4, 5)]
  rownames(output) <- substr(rownames(output), 2, nchar(rownames(output)))

  print(output)
  cat("\n\n")
}
#' @export
print.hJAM_gprior <- function(x, ...) {
  cat("--------------------------------------------", "\n")
  cat("           hJAM gprior output              ", "\n")
  cat("--------------------------------------------", "\n")

  cat(paste0("Number of SNPs used in model: ", x$numSNP), "\n\n")

  output =  as.data.frame(unclass(x))[c(3, 4, 5)]
  rownames(output) <- substr(rownames(output), 2, nchar(rownames(output)))

  print(output)
  cat("\n\n")
}
#' @export
print.hJAM_egger <- function(x, ...) {
  cat("--------------------------------------------", "\n")
  cat("             hJAM egger output              ", "\n")
  cat("--------------------------------------------", "\n")

  cat(paste0("Number of SNPs used in model: ", x$numSNP), "\n\n")

  cat("Exposures\n")
  exp_output = as.data.frame(unclass(x))[c(3, 4, 5)]
  rownames(exp_output) <- substr(rownames(exp_output), 2, nchar(rownames(exp_output)))
  print(exp_output)

  cat("\nIntercept\n")
  int_output = cbind(x$Est.Int, x$StdErr.Int, x$Pvalue.Int)
  colnames(int_output) = c("Est.Int", "StdErr.Int", "Pvalue.Int")
  print(int_output)
  cat("\n\n")
}
