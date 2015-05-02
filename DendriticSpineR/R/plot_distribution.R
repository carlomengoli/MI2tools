#' Function plots crossed effects
#'
#' Function plots densities for given varialbe divided into groups (f1) plotted on different panels (strat)
#'
#' @param data data.frame with data
#' @param var variable of interest
#' @param f1 first effect (used as fill)
#' @param strat stratification (Animal)
#'
#' @import ggplot2
#'
#' @export

plot_distribution <- function(data, var, f1="Group", strat = "Animal") {
  quant <- quantile(data[,var], c(0.001, 0.999))

  ggplot(data, aes_string(fill=f1, x=var)) +
    geom_density(adjust=1, alpha=0.5) +
    coord_cartesian(xlim=quant) +
    facet_wrap(as.formula(paste0("~", strat)))
}
