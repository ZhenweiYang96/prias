library(ggplot2)
library(ggpubr)

source("src/clinical_gap3/plots/npmleCurve.R")
load("Rdata/gap3/PRIAS_2019/validation/auc_df.Rdata")

FONT_SIZE = 14

auc_recalibdf = auc_df[auc_df$model=="Recalibrated",]
auc_recalibdf = auc_recalibdf[!auc_recalibdf$cohort %in% "UCSF" | auc_recalibdf$t_horiz<9,]
auc_recalibplot = ggplot(data=auc_recalibdf) + 
  geom_line(aes(x=t_horiz, y=mean, group=cohort, color=cohort)) + 
  scale_color_manual(values = colormap)+
  scale_x_continuous(breaks = seq(1,9, by=1), limits = c(1,9.1)) + 
  theme_bw() + theme(text=element_text(size=FONT_SIZE), legend.title = element_blank()) +
  ylab("AUC (higher is better)") + xlab("Follow-up time (years)") +
  ylim(0.5,1)

#calibration plot
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
                      cum_risk=c(calib_Hopkins, calib_KCL,
                                 calib_MSKCC, calib_Toronto, 
                                 calib_MUSIC, calib_PRIAS, calib_UCSF),
                      cohort=rep(c("Hopkins", "KCL", "MSKCC", "Toronto", 
                                   "MUSIC", "PRIAS", "UCSF"), each=length(calib_pred_times)))

calib_df$max_time = sapply(calib_df$cohort, function(x){
  reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==x]
})

calib_df = calib_df[calib_df$pred_time <= calib_df$max_time,]

before_recalib_main_manuscript = npmle_plot_all +
  geom_line(data=calib_df, aes(x=pred_time, y=cum_risk, color=cohort), 
            linetype="dashed")+
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohortnames,
                 fill=cohortnames), color='white')

auc_calib_plot = ggarrange(auc_recalibplot, before_recalib_main_manuscript, nrow=1, ncol = 2,
                           align = "h", common.legend = T, legend = "bottom", labels = "AUTO")
ggsave(auc_calib_plot, device = cairo_ps,
       file="report/clinical/images/auc_beforecalib.eps", width = 8, height=5.5)

