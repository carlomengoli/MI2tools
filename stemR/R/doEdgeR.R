#' Function performs edgeR two-sample alculations
#'
#' @param counts numeric matrix of read counts.
#' @param group vector or factor giving the experimental group/condition for each sample/library.
#'
#' @import edgeR
#'
#' @export

doEdgeR <- function(counts, group) {
  y <- DGEList(counts=counts, group=group)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  et <- exactTest(y)
  list(counts = y,
       tests = et)
}
