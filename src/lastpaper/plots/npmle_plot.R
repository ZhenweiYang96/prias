library(JMbayes)
load("Rdata/gap3/PRIAS_2019/npmle_all.Rdata")

npmle_plotdf_all=do.call('rbind', 
                         lapply(c("PRIAS"), FUN = function(name){
                           survProb = 1 - cumsum(npmle_all[[name]]$pf)
                           survProb = c(1, survProb)
                           
                           survIntervals = npmle_all[[name]]$intmap
                           survIntervals = cbind(c(0,0), survIntervals)
                           
                           timePoints = as.numeric(survIntervals)
                           survProbs = c(1,as.numeric(rep(survProb, each=2)))[1:length(timePoints)]
                           
                           return(data.frame('Cohort'=name,timePoints=timePoints, riskProbs=1-survProbs))
                         }))


FONT_SIZE=14
npmle_plot = ggplot() + 
  geom_line(aes(x=npmle_plotdf_all$timePoints, 
                y=npmle_plotdf_all$riskProbs)) +  
  coord_cartesian(xlim=c(0,10)) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        legend.position = "none",
        legend.text = element_text(size=FONT_SIZE-4),
        axis.line = element_line())+
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"),
                     limits = c(0,1)) + 
  ylab("Cumulative-risk of progression (%)") +
  xlab("Follow-up time (years)")

ggsave(npmle_plot, filename = "report/lastpaper/images/npmle_plot.eps",
       device = cairo_ps,  height = 5.5, width=5.5)
