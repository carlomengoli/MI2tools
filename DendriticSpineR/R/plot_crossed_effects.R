#' Function plots crossed effects
#'
#' Function creates folder with downloaded subtitles (as .zip file or unzipped .txt or .srt file)
#'
#' @param data data.frame with data
#' @param var variable of interest
#' @param trans transformation of the var variable before ANOVA
#' @param inv inverse transformation of the trans(var)
#' @param f1 first effect
#' @param f2 second effect
#' @param strat stratification (Animal)
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export

plot_crossed_effects <- function(data, var, trans = I, inv = I, f1="group", f2="condition", strat = "Animal") {
  ndata <- data
  ndata[,var] <- trans(data[,var])
  form <- as.formula(paste0(var, " ~ ", f1, ":", f2, "/factor(",strat,")-1"))
  model <- summary(lm(form, ndata))

  co2 <- data.frame(n = rownames(model$coef), model$coef[,1:2])
  co2 <- co2[!grepl(co2$n, pattern = "factor(", fixed = TRUE),]
  co2 <- co2 %>% dplyr::mutate(l1 = inv(Estimate - Std..Error),
                  l2 = inv(Estimate),
                  l3 = inv(Estimate + Std..Error)) %>%
    extract(n, c(f1, f2), paste0(f1, "([[:alnum:]]+):", f2, "([[:alnum:]]+)"))

  ggplot(co2, aes_string(x=f1, color=f2, ymin="l1", y="l2", ymax="l3")) +
    geom_errorbar(width=0.1, position=position_dodge(width=0.2)) +
    geom_point(size=5, position=position_dodge(width=0.2)) +
    theme_bw() + ylab(var)
}
