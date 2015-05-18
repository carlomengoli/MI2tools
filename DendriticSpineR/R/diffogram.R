diffogram <- function(lsmodel) {
  # effects
  tmp1 <- confint(lsmodel$lsmeans)
  effects <- data.frame(labels = tmp1[,attr(tmp1, "pri.vars")],
                   values = tmp1[,attr(tmp1, "estName")])
  rownames(effects) <- effects$labels
  
  # confidence intervals
  tmp2 <- confint(ms6$contrasts)
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
                        wsp_d1 = ci[,"values"], 
                        wsp_d2 = ci[,"ciL"], 
                        wsp_d3 = ci[,"ciR"])
  
    
}




ggplot(to_plot, aes(x=wsp_x, y=wsp_y)) + 
  geom_hline(data=ef, aes(yintercept=values), lty=2, color="grey") + 
  geom_vline(data=ef, aes(xintercept=values), lty=2, color="grey") + 
  geom_point(size=2) + 
  geom_segment(aes(x=wsp_x + wsp_d2, xend=wsp_x + wsp_d3,
                   y=wsp_y - wsp_d2, yend=wsp_y - wsp_d3)) + 
  geom_abline(intercept=0, slope=1) 


