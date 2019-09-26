library(JMbayes)
load("Rdata/gap3/PRIAS_2019/npmle_all.Rdata")

npmle_plotdf_all=do.call('rbind', 
                         lapply(c("Hopkins", "London-KCL", "MSKCC", "PRIAS", "Toronto", "MUSIC"), FUN = function(name){
                           survProb = 1 - cumsum(npmle_all[[name]]$pf)
                           survProb = c(1, survProb)
                           
                           survIntervals = npmle_all[[name]]$intmap
                           survIntervals = cbind(c(0,0), survIntervals)
                           
                           timePoints = as.numeric(survIntervals)
                           survProbs = c(1,as.numeric(rep(survProb, each=2)))[1:length(timePoints)]
                           
                           return(data.frame('Cohort'=name,timePoints=timePoints, riskProbs=1-survProbs))
                         }))
levels(npmle_plotdf_all$Cohort)[2] = c("KCL")

npmle_plotdf_all$time_10pat_risk_set = sapply(npmle_plotdf_all$Cohort, function(x){
  reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==x]
})
npmle_plotdf_all = npmle_plotdf_all[npmle_plotdf_all$timePoints <= npmle_plotdf_all$time_10pat_risk_set,]

cohort_names = unique(npmle_plotdf_all$Cohort)
cohort_labpos_x = as.numeric(by(npmle_plotdf_all$Cohort, data = npmle_plotdf_all$timePoints, max))
cohort_labpos_y = as.numeric(by(npmle_plotdf_all$Cohort, data = npmle_plotdf_all$riskProbs, max))

FONT_SIZE=13
npmle_plot_all = ggplot() + 
  geom_line(aes(x=npmle_plotdf_all$timePoints, 
                y=npmle_plotdf_all$riskProbs, 
                group=npmle_plotdf_all$Cohort, 
                color=npmle_plotdf_all$Cohort)) +  
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohort_names,
                 fill=cohort_names), color='white')+
  coord_cartesian(xlim=c(0,8)) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        axis.text=element_text(size=FONT_SIZE),
        legend.position = "none",
        legend.text = element_text(size=FONT_SIZE-4),
        axis.line = element_line(), 
        plot.margin = margin(0, 0, 0, 0, "pt")) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"),
                     limits = c(0,1)) + 
  ylab("Cumulative risk of reclassification (%)") +
  xlab("Follow-up time (years)")

ggsave(filename = "report/clinical/images/npmle_plot.eps",
       plot=npmle_plot_all, device=cairo_ps, height=5.5, width=6, dpi = 500)
