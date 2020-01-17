library(ggplot2)

load("Rdata/gap3/PRIAS_2019/motherdata.Rdata")

load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/Hopkins.Rdata")
diff_Hopkins= c(cumrisk_models[[2]] - cumrisk_models[[1]])
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/Toronto.Rdata")
diff_Toronto= c(cumrisk_models[[3]] - cumrisk_models[[1]])
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/MSKCC.Rdata")
diff_MSKCC= c(cumrisk_models[[3]] - cumrisk_models[[1]])
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/KCL.Rdata")
diff_KCL= c(cumrisk_models[[3]] - cumrisk_models[[1]])
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/MUSIC.Rdata")
diff_MUSIC= c(cumrisk_models[[3]] - cumrisk_models[[1]])
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/UCSF.Rdata")
diff_UCSF= c(cumrisk_models[[3]] - cumrisk_models[[1]])

calib_pred_times = seq(0, 10, 0.1)
diff_df = data.frame(cohort = rep(cohortnames[cohortnames!="PRIAS"], sapply(list(diff_Hopkins, diff_Toronto, diff_MSKCC, 
                                                           diff_MUSIC, diff_KCL, diff_UCSF), length)),
                     time = calib_pred_times,
                     diff = c(diff_Hopkins, diff_Toronto, diff_MSKCC, diff_MUSIC, diff_KCL, diff_UCSF))

diff_df$max_time = sapply(diff_df$cohort, function(x){
  reclassification_df$time_10pat_risk_set[as.character(reclassification_df$Cohort)==as.character(x)]
})
diff_df = diff_df[diff_df$time <= diff_df$max_time,]

diff_plot = ggplot(diff_df) + geom_boxplot(aes(y=diff,x=cohort), outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-0.5,0.5,by=0.1), 
                     limits = c(-0.5,0.5),
                     labels = paste0(seq(-0.5,0.5,by = 0.1)*100, "%")) + theme_bw() +
  theme(text=element_text(size=18)) +
  xlab("Cohort") + ylab("Difference in predicted\ncumulative-risk")

ggsave(diff_plot, device = cairo_ps,
       file="report/clinical/images/calib_insmall_after.eps")
