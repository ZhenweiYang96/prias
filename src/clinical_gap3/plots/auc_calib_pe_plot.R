library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/validation/auc_df.Rdata")
load("Rdata/gap3/PRIAS_2019/validation/pe_df.Rdata")

FONT_SIZE = 14

auc_origdf = auc_df[auc_df$model=="Original",]
auc_recalibdf = auc_df[auc_df$model=="Recalibrated",]

auc_origplot = ggplot(data=auc_origdf) + 
  geom_line(aes(x=t_horiz, y=mean, group=cohort, color=cohort)) + 
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1)) + 
  scale_color_manual(values = colormap)+
  geom_hline(yintercept = 0.5, linetype='dashed') + 
  geom_label(aes(x=4.5, y=0.4), label="0.5 is the AUC\nfor random discrimination", 
             size=3.5) +
  scale_x_continuous(breaks = seq(1,8, by=1), limits = c(1,8)) + 
  theme_bw() + theme(text=element_text(size=FONT_SIZE), legend.title = element_blank()) +
  ylab("AUC (higher is better)") + xlab("Follow-up time (years)")

auc_recalibplot = ggplot(data=auc_recalibdf) + 
  geom_line(aes(x=t_horiz, y=mean, group=cohort, color=cohort)) + 
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1)) + 
  scale_color_manual(values = colormap)+
  geom_hline(yintercept = 0.5, linetype='dashed') + 
  geom_label(aes(x=4.5, y=0.4), label="0.5 is the AUC\nfor random discrimination", 
             size=3.5) +
  scale_x_continuous(breaks = seq(1,8, by=1), limits = c(1,8)) + 
  theme_bw() + theme(text=element_text(size=FONT_SIZE), legend.title = element_blank()) +
  ylab("AUC (higher is better)") + xlab("Follow-up time (years)")

mapeplot = ggplot(data=pe_df) + 
  geom_line(aes(x=t_horiz, y=mean_mape, group=cohort, color=cohort)) + 
  scale_color_manual(values = colormap)+
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1)) + 
  scale_x_continuous(breaks = seq(1,8, by=1), limits = c(1,8)) + 
  theme_bw() + theme(text=element_text(size=FONT_SIZE), legend.title = element_blank()) +
  ylab("MAPE (lower is better)") + xlab("Follow-up time (years)")

auc_pe = ggarrange(auc_recalibplot, mapeplot, nrow=1, ncol = 2, 
                   common.legend = T, legend = "bottom",labels = "AUTO")

ggsave(auc_pe, device = cairo_ps,
       file="report/clinical/images/auc_pe_recalib.eps", width = 7, height=5.5)


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

calib_df = data.frame(pred_time=calib_pred_times,
                      cum_risk=c(calib_Hopkins, calib_KCL,
                                 calib_MSKCC, calib_Toronto, calib_MUSIC, calib_PRIAS),
                      cohort=rep(c("Hopkins", "KCL", "MSKCC", "Toronto", "MUSIC", "PRIAS"), each=length(calib_pred_times)))

calib_df$max_time = sapply(calib_df$cohort, function(x){
  reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==x]
})

calib_df = calib_df[calib_df$pred_time <= calib_df$max_time,]

before_recalib = npmle_plot_all +
  geom_line(data=calib_df, aes(x=pred_time, y=cum_risk, color=cohort), 
            linetype="dashed")+
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohort_names,
                 fill=cohort_names), color='white') +
  ggtitle("Before recalibration")

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
load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/PRIAS.Rdata")
calib_PRIAS= rowMeans(cumrisk_models[[2]],na.rm = T)

recalib_df = data.frame(pred_time=calib_pred_times,
                      cum_risk=c(calib_Hopkins, calib_KCL,
                                 calib_MSKCC, calib_Toronto, calib_MUSIC, calib_PRIAS),
                      cohort=rep(c("Hopkins", "KCL", "MSKCC", "Toronto", "MUSIC", "PRIAS"), each=length(calib_pred_times)))

recalib_df$max_time = sapply(recalib_df$cohort, function(x){
  reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==x]
})

recalib_df = recalib_df[recalib_df$pred_time <= recalib_df$max_time,]

after_recalib = npmle_plot_all +
  geom_line(data=recalib_df, aes(x=pred_time, y=cum_risk, color=cohort), 
            linetype="dashed")+
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohort_names,
                 fill=cohort_names), color='white') +
  ggtitle("   After recalibration") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

calib_in_large =  ggarrange(before_recalib, after_recalib, nrow=1, ncol = 2,
                             align = "h", widths = c(1.25,1),
                             common.legend = F, legend = "none", labels = "AUTO")

ggsave(calib_in_large, device = cairo_ps,
       file="report/clinical/images/calib_before_after.eps", width = 8, height=5.5)


#now auc calib plot
after_recalib_main_manuscript = npmle_plot_all +
  geom_line(data=recalib_df, aes(x=pred_time, y=cum_risk, color=cohort), 
            linetype="dashed")+
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohort_names,
                 fill=cohort_names), color='white') +
  theme(text=element_text(size=FONT_SIZE))
auc_aftercalib_plot = ggarrange(auc_recalibplot, after_recalib_main_manuscript, nrow=1, ncol = 2,
                           align = "h", common.legend = T, legend = "bottom", labels = "AUTO")
ggsave(auc_aftercalib_plot, device = cairo_ps,
       file="report/clinical/images/auc_aftercalib.eps", width = 8, height=5.5)

before_recalib_main_manuscript = npmle_plot_all +
  geom_line(data=calib_df, aes(x=pred_time, y=cum_risk, color=cohort), 
            linetype="dashed")+
  geom_label(aes(x=cohort_labpos_x, 
                 y=cohort_labpos_y, 
                 label=cohort_names,
                 fill=cohort_names), color='white') +
  theme(text=element_text(size=FONT_SIZE))
auc_beforecalib_plot = ggarrange(auc_origplot, before_recalib_main_manuscript, nrow=1, ncol = 2,
                           align = "h", common.legend = T, legend = "bottom", labels = "AUTO")
ggsave(auc_beforecalib_plot, device = cairo_ps,
       file="report/clinical/images/auc_beforecalib.eps", width = 8, height=5.5)


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

calib_pred_times = seq(0, 10, 0.1)
diff_df = data.frame(cohort = rep(c("Hopkins", "Toronto", "MSKCC", "KCL", "MUSIC"), 
                                  sapply(list(diff_Hopkins, diff_Toronto, diff_MSKCC, diff_KCL, diff_MUSIC), length)),
                     time = calib_pred_times,
  diff = c(diff_Hopkins, diff_Toronto, diff_MSKCC, diff_KCL, diff_MUSIC))

diff_df$max_time = sapply(diff_df$cohort, function(x){
  reclassification_df$time_10pat_risk_set[as.character(reclassification_df$Cohort)==as.character(x)]
})
diff_df = diff_df[diff_df$time <= diff_df$max_time,]

diff_plot = ggplot(diff_df) + geom_boxplot(aes(y=diff,x=cohort), outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,1,by=0.25), limits = c(-1,1),
                     labels = paste0(seq(-1,1,by = 0.25)*100, "%")) + theme_bw() +
  theme(text=element_text(size=15)) +
  xlab("Cohort") + ylab("Difference in predicted cumulative-risk")

ggsave(diff_plot, device = cairo_ps,
       file="report/clinical/images/calib_insmall_after.eps")
