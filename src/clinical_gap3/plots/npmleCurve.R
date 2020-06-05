library(JMbayes)
library(icenReg)
load("Rdata/gap3/PRIAS_2019/motherdata.Rdata")

getSurvProbNPMLE = function(pred.time, baseline_est_int){
  #turnbull intervals are hard to understand
  index2 = tail(which(baseline_est_int[,2]<=pred.time),1)
  if(baseline_est_int[index2 + 1, 1]<=pred.time){
    pred.time.surv = baseline_est_int[index2 + 1,3]
  }else{
    pred.time.surv = baseline_est_int[index2, 3]
  }
  
  return(pred.time.surv)
}

npmle_all = vector("list", length = length(cohortnames))
for(cohort in cohortnames){
  longds.id = get(paste0("longdata_", cohort, ".id"))
  npmle_all[[cohort]] = ic_sp(formula = Surv(time = latest_survival_time, time2 = earliest_failure_time,
                                             type = "interval2")~1, data = longds.id)
}

npmle_plotdf_all=do.call('rbind', lapply(cohortnames, FUN = function(cohort){
  
  baseline_est = getSCurves(npmle_all[[cohort]])
  baseline_est_int = baseline_est$Tbull_ints
  baseline_est_int = cbind(baseline_est_int, baseline_est$S_curves$baseline)
  baseline_est_int = rbind(c(0,0,1), baseline_est_int)
  
  max_time = reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==cohort]
  timePoints = seq(0, max_time, length.out = 100)
  survProbs = sapply(timePoints, FUN = getSurvProbNPMLE, baseline_est_int=baseline_est_int)
  
  return(data.frame('Cohort'=cohort, timePoints=timePoints, riskProbs=1-survProbs))
}))


npmle_plotdf_all_homogenizationclaim=do.call('rbind', lapply(cohortnames, FUN = function(cohort){
  
  baseline_est = getSCurves(npmle_all[[cohort]])
  baseline_est_int = baseline_est$Tbull_ints
  baseline_est_int = cbind(baseline_est_int, baseline_est$S_curves$baseline)
  baseline_est_int = rbind(c(0,0,1), baseline_est_int)
  
  max_time = reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==cohort]
  timePoints = seq(1, max_time, length.out = 100)
  survProbs = sapply(timePoints, FUN = getSurvProbNPMLE, baseline_est_int=baseline_est_int)
  
  survProb_yearone = getSurvProbNPMLE(pred.time = 1, baseline_est_int)
  survProbs = survProbs/survProb_yearone
  
  return(data.frame('Cohort'=cohort, timePoints=timePoints, riskProbs=1-survProbs))
}))

cohort_labpos_x = as.numeric(by(npmle_plotdf_all$Cohort, data = npmle_plotdf_all$timePoints, max))
cohort_labpos_y = as.numeric(by(npmle_plotdf_all$Cohort, data = npmle_plotdf_all$riskProbs, max))

FONT_SIZE=14
npmle_plot_all = ggplot() + 
  geom_line(aes(x=npmle_plotdf_all$timePoints, 
                y=npmle_plotdf_all$riskProbs, 
                group=npmle_plotdf_all$Cohort, 
                color=npmle_plotdf_all$Cohort)) +  
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohortnames,
                 fill=cohortnames), color='white')+
  scale_color_manual(values=colormap)+
  scale_fill_manual(values=colormap)+
  scale_x_continuous(breaks=0:9, limits = c(0,9)) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        legend.position = "none",
        legend.text = element_text(size=FONT_SIZE-4),
        axis.line = element_line())+
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"),
                     limits = c(0,1)) + 
  ylab("Cause-specific cumulative upgrading-risk (%)") +
  xlab("Follow-up time (years)")
print(npmle_plot_all)

ggsave(filename = "report/clinical/images/npmle_plot.eps",
       plot=npmle_plot_all, device=cairo_ps, height=5.5, width=6, dpi = 500)


