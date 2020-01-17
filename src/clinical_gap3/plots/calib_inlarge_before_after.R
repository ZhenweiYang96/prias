library(ggplot2)
library(ggpubr)

FONT_SIZE = 14

#calibration plot
source("src/clinical_gap3/plots/npmleCurve.R")
calib_pred_times = seq(0, 10, 0.1)

load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/Hopkins.Rdata")
calib_Hopkins= rowMeans(cumrisk_models[[2]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/Toronto.Rdata")
calib_Toronto= rowMeans(cumrisk_models[[2]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/MSKCC.Rdata")
calib_MSKCC= rowMeans(cumrisk_models[[2]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/KCL.Rdata")
calib_KCL= rowMeans(cumrisk_models[[2]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/MUSIC.Rdata")
calib_MUSIC= rowMeans(cumrisk_models[[2]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/PRIAS.Rdata")
calib_PRIAS= rowMeans(cumrisk_models[[2]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/UCSF.Rdata")
calib_UCSF= rowMeans(cumrisk_models[[2]],na.rm = T)

calib_df = data.frame(pred_time=calib_pred_times,
                      cum_risk=c(calib_Hopkins, calib_Toronto,
                                 calib_MSKCC, calib_MUSIC,
                                 calib_KCL, calib_PRIAS, calib_UCSF),
                      cohort=rep(cohortnames, each=length(calib_pred_times)))

calib_df$max_time = sapply(calib_df$cohort, function(x){
  reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==x]
})

calib_df = calib_df[calib_df$pred_time <= calib_df$max_time,]

before_recalib = npmle_plot_all +
  geom_line(data=calib_df, aes(x=pred_time, y=cum_risk, color=cohort), 
            linetype="dashed")+
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohortnames,
                 fill=cohortnames), color='white') +
  ggtitle("Before recalibration")

####################
# After
####################

load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/PRIAS.Rdata")
calib_PRIAS= rowMeans(cumrisk_models[[2]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/Hopkins.Rdata")
calib_Hopkins= rowMeans(cumrisk_models[[2]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/Toronto.Rdata")
calib_Toronto= rowMeans(cumrisk_models[[3]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/MSKCC.Rdata")
calib_MSKCC= rowMeans(cumrisk_models[[3]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/KCL.Rdata")
calib_KCL= rowMeans(cumrisk_models[[3]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/MUSIC.Rdata")
calib_MUSIC= rowMeans(cumrisk_models[[3]],na.rm = T)
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/UCSF.Rdata")
calib_UCSF= rowMeans(cumrisk_models[[3]],na.rm = T)

recalib_df = data.frame(pred_time=calib_pred_times,
                      cum_risk=c(calib_Hopkins, calib_Toronto,
                                 calib_MSKCC, calib_MUSIC,
                                 calib_KCL, calib_PRIAS, calib_UCSF),
                      cohort=rep(cohortnames, each=length(calib_pred_times)))

recalib_df$max_time = sapply(recalib_df$cohort, function(x){
  reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==x]
})

recalib_df = recalib_df[recalib_df$pred_time <= recalib_df$max_time,]

after_recalib = npmle_plot_all +
  geom_line(data=recalib_df, aes(x=pred_time, y=cum_risk, color=cohort), 
            linetype="dashed")+
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohortnames,
                 fill=cohortnames), color='white') +
  ggtitle("   After recalibration") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

calib_in_large =  ggarrange(before_recalib, after_recalib, nrow=1, ncol = 2,
                            align = "h", widths = c(1.25,1),
                            common.legend = F, legend = "none", labels = "AUTO")

ggsave(calib_in_large, device = cairo_ps,
       file="report/clinical/images/calib_before_after.eps", width = 8, height=5.5)
