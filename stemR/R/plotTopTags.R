#' Function plots top genes
#'
#' @param counts matrix (small) for top genes.
#' @param groups groups for following collumns in counts.
#'
#' @import ggplot2
#' @import tidyr
#'
#' @export

plotTopTags <- function(counts, groups) {
  topT <- as.data.frame(counts)
  vn <- make.names(groups, unique = TRUE)
  names(groups) <- vn
  colnames(topT) <- vn
  topT$gene <- rownames(topT)
  topTlong <- gather(topT, group, value, -gene)
  topTlong$group <- groups[as.character(topTlong$group)]

  ggplot(topTlong, aes(x=group, y=value, color=group)) +
    geom_boxplot() +
    geom_point(position=position_jitter(0.2,0)) +
    facet_wrap(~gene, scales = "free_y") + theme_bw()
}
