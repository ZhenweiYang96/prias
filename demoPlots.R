npmle_plotdf_all=do.call('rbind', 
                         lapply(c("Hopkins", "London-KCL", "MSKCC", "PRIAS", "Toronto"), FUN = function(name){
                           survProb = 1 - cumsum(npmle_all[[name]]$pf)
                           survProb = c(1, survProb)
                           
                           survIntervals = npmle_all[[name]]$intmap
                           survIntervals = cbind(c(0,0), survIntervals)
                           
                           timePoints = as.numeric(survIntervals)
                           survProbs = c(1,as.numeric(rep(survProb, each=2)))[1:length(timePoints)]
                           
                           return(data.frame('Cohort'=name,timePoints=timePoints, riskProbs=1-survProbs))
                         }))
levels(npmle_plotdf_all$Cohort)[1:2] = c("JHAS", "KCL")

FONT_SIZE=13
npmle_plot_all = ggplot(data=npmle_plotdf_all) + 
  geom_line(aes(x=timePoints, y=riskProbs, group=Cohort, 
                color=Cohort)) +  
  coord_cartesian(xlim=c(0,10)) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        axis.text=element_text(size=FONT_SIZE),
        legend.position = "bottom",
        legend.text = element_text(size=FONT_SIZE-4),
        axis.line = element_line(), 
        plot.margin = margin(0, 0, 0, 0, "pt")) + 
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"),
                     limits = c(0,1)) + 
  ylab("Cumulative risk of reclassification (%)") +
  xlab("Follow-up time (years)")

ggsave(filename = "report/clinical/images/npmle_plot.eps",
       plot=npmle_plot_all, device=cairo_ps, height=5.5, width=5.5, dpi = 500)

