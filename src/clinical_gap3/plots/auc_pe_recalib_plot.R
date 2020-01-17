library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/validation/auc_df.Rdata")
load("Rdata/gap3/PRIAS_2019/validation/pe_df.Rdata")

FONT_SIZE = 14
# auc_origdf = auc_df[auc_df$model=="Original",]
# auc_origplot = ggplot(data=auc_origdf) + 
#   geom_line(aes(x=t_horiz, y=mean, group=cohort, color=cohort)) + 
#   scale_color_manual(values = colormap)+
#   scale_x_continuous(breaks = seq(1,8, by=1), limits = c(1,8.1)) + 
#   theme_bw() + theme(text=element_text(size=FONT_SIZE), legend.title = element_blank()) +
#   ylab("AUC (higher is better)") + xlab("Follow-up time (years)") +
#   ylim(0.5,1)

auc_recalibdf = auc_df[auc_df$model=="Recalibrated",]
auc_recalibdf = auc_recalibdf[!auc_recalibdf$cohort %in% "UCSF" | auc_recalibdf$t_horiz<9,]
auc_recalibplot = ggplot(data=auc_recalibdf) + 
  geom_line(aes(x=t_horiz, y=mean, group=cohort, color=cohort)) + 
  scale_color_manual(values = colormap)+
  scale_x_continuous(breaks = seq(1,9, by=1), limits = c(1,9.1)) + 
  theme_bw() + theme(text=element_text(size=FONT_SIZE), legend.title = element_blank()) +
  ylab("AUC (higher is better)") + xlab("Follow-up time (years)") +
  ylim(0.5,1)

pe_df = pe_df[!pe_df$cohort %in% "UCSF" | pe_df$t_horiz<9,]
mapeplot = ggplot(data=pe_df) + 
  geom_line(aes(x=t_horiz, y=mean_mape, group=cohort, color=cohort)) + 
  scale_color_manual(values = colormap)+
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1)) + 
  scale_x_continuous(breaks = seq(1,9, by=1), limits = c(1,9.1)) + 
  theme_bw() + theme(text=element_text(size=FONT_SIZE), legend.title = element_blank()) +
  ylab("MAPE (lower is better)") + xlab("Follow-up time (years)")

auc_pe = ggarrange(auc_recalibplot, mapeplot, nrow=1, ncol = 2, 
                   common.legend = T, legend = "bottom",labels = "AUTO")

ggsave(auc_pe, device = cairo_ps,
       file="report/clinical/images/auc_pe_recalib.eps", width = 7, height=5.5)