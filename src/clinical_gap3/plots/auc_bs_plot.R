library(ggplot2)
library(ggpubr)

files = list.files("Rdata/gap3/PRIAS_2019/auc_prederr/", full.names = T)

t_horizs = seq(1, 10, 0.5)
bs_number = 1:25

cohortNames = c("JHAS", "KCL", "MSKCC", "MUSIC",  "PRIAS", "Toronto")

auc_rmspe_df = do.call('rbind',lapply(1:length(files), function(i){
  load(files[[i]])
  
  aucs = as.numeric(sapply(auc_prederr_bs, function(bs){
    sapply(bs$auc_list, "[[", "auc")
  }))
  
  prederrs = as.numeric(sapply(auc_prederr_bs, function(bs){
    sapply(bs$prederr_list, "[[", "prederr")
  }))
  
  return(data.frame(Cohort=cohortNames[i],
                    bs_number=rep(bs_number, each=length(t_horizs)),
                    t_start = t_horizs-1, t_horiz=t_horizs,
                    auc=aucs, rmspe=sqrt(prederrs)))
}))

plotDf=do.call('rbind', by(auc_rmspe_df, INDICES = auc_rmspe_df$Cohort, FUN = function(cohort_auc_rmspe_df){
  mean_auc = c(by(cohort_auc_rmspe_df$auc, cohort_auc_rmspe_df$t_start, mean, na.rm=T))
  lower_auc = c(by(cohort_auc_rmspe_df$auc, cohort_auc_rmspe_df$t_start, quantile, probs=0.025, na.rm=T))
  upper_auc = c(by(cohort_auc_rmspe_df$auc, cohort_auc_rmspe_df$t_start, quantile, probs=0.975, na.rm=T))

  mean_rmspe = c(by(cohort_auc_rmspe_df$rmspe, cohort_auc_rmspe_df$t_start, mean, na.rm=T))
  lower_rmspe = c(by(cohort_auc_rmspe_df$rmspe, cohort_auc_rmspe_df$t_start, quantile, probs=0.025, na.rm=T))
  upper_rmspe = c(by(cohort_auc_rmspe_df$rmspe, cohort_auc_rmspe_df$t_start, quantile, probs=0.975, na.rm=T))
  
  return(data.frame(Cohort = cohort_auc_rmspe_df$Cohort[1],
                    t_start = t_horizs-1, t_horiz=t_horizs,
                    mean_auc,lower_auc,upper_auc,
                    mean_rmspe,lower_rmspe,upper_rmspe))
}))

FONT_SIZE = 15

aucplot = ggplot(data=plotDf) + 
  geom_line(aes(x=t_horiz, y=mean_auc, group=Cohort, color=Cohort)) + 
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1)) + 
  geom_hline(yintercept = 0.5, linetype='dashed') + 
  geom_label(aes(x=3, y=0.45), label="0.5 is the AUC for random discrimination", 
             size=3.5) +
  scale_x_continuous(breaks = seq(1,5, by=1), limits = c(1,5)) + 
  theme_bw() + theme(text=element_text(size=FONT_SIZE)) +
  ylab("AUC (higher is better)") + xlab("Follow-up time (years)")

rmspeplot = ggplot(data=plotDf) + 
  geom_line(aes(x=t_horiz, y=mean_rmspe, group=Cohort, color=Cohort)) + 
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1)) + 
  scale_x_continuous(breaks = seq(1,5, by=1), limits = c(1,5)) + 
  theme_bw() + theme(text=element_text(size=FONT_SIZE)) +
  ylab("RMSPE (lower is better)") + xlab("Follow-up time (years)")

ggsave(ggarrange(aucplot, rmspeplot, nrow=1, ncol = 2,
          common.legend = T, legend = "bottom", labels = "AUTO"), device = cairo_ps,
       file="report/clinical/images/auc_pe.eps", width = 8, height=5.5)
