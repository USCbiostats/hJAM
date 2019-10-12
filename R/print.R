#' @export
print.hJAM <- function(x, ...) {
  cat("hJAM output", "\n")
  cat("-----------", "\n\n")

  matrix <- as.data.frame(unclass(x))
  rownames(matrix) <- paste0("Gx ", rownames(matrix), ":")
  print(matrix)
}
