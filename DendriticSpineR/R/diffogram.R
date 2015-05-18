diffogram <- function(lsmodel) {
  # effects
  tmp1 <- confint(lsmodel$lsmeans)
  effects <- data.frame(labels = tmp1[,attr(tmp1, "pri.vars")],
                   values = tmp1[,attr(tmp1, "estName")])
  rownames(effects) <- effects$labels
  
  # confidence intervals
  tmp2 <- confint(lsmodel$contrasts)
  ci <- data.frame(labels = tmp2[,attr(tmp2, "pri.vars")],
                    values = tmp2[,attr(tmp2, "estName")],
                    ciL = tmp2[,attr(tmp2, "clNames")[1]],
                    ciR = tmp2[,attr(tmp2, "clNames")[2]])
  
  # indexes
  ll <- strsplit(as.character(ci$labels), split=" - ")
  ci_names <- data.frame(
    ci_left_name  = sapply(ll, `[`, 1),
    ci_right_name = sapply(ll, `[`, 2),
    stringsAsFactors = FALSE
  )

  to_plot <- data.frame(name_x = ci_names$ci_left_name,
                        name_y = ci_names$ci_right_name,
                        wsp_x = effects[ci_names$ci_left_name,"values"], 
                        wsp_y = effects[ci_names$ci_right_name,"values"],
                        wsp_x_y = ci[,"values"], 
                        wsp_x_y_ci_left = ci[,"ciL"] - ci[,"values"], 
                        wsp_x_y_ci_right = ci[,"ciR"] - ci[,"values"],
                        significant = ifelse(ci[,"ciL"] * ci[,"ciR"] > 0,
                                             "significant", "non-significant"))

  # ranges
  spec <- range(effects$values) + max(c(abs(to_plot$wsp_x_y_ci_left), abs(to_plot$wsp_x_y_ci_right))) * c(-0.5,0.5)
  
  # the plot
  ggplot(to_plot, aes(x=wsp_x, y=wsp_y)) + 
    geom_hline(data=effects, aes(yintercept=values), lty=2, color="grey") + 
    geom_vline(data=effects, aes(xintercept=values), lty=2, color="grey") + 
    geom_point(size=2) + 
    geom_segment(aes(x=wsp_y + wsp_x_y - wsp_x_y_ci_left/2, 
                     xend=wsp_y + wsp_x_y - wsp_x_y_ci_right/2,
                     y=wsp_y  + wsp_x_y_ci_left/2, 
                     yend=wsp_y + wsp_x_y_ci_right/2,
                     color=significant,
                     lty=significant)) + 
    geom_abline(intercept=0, slope=1) +
    xlim(spec) + ylim(spec) +
    scale_color_manual(values=c("navyblue", "red4")) +
    scale_linetype_manual(values=c(2,1)) +
    theme(panel.background	= element_rect(fill = "white"))
  
    
}






