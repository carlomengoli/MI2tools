#' Function performs filtering, removes rows with low counts
#'
#' @param counts numeric matrix of read counts.
#' @param avg.min minimum average per row.
#'
#' @export

doFiltering <- function(counts, avg.min = 1) {
  counts[rowMeans(counts) >= avg.min,,drop=FALSE]
}
