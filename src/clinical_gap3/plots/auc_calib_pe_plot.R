library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/validation/auc_df.Rdata")
load("Rdata/gap3/PRIAS_2019/validation/pe_df.Rdata")

FONT_SIZE = 15

auc_df = auc_df[auc_df$model=="Recalibrated",]

aucplot = ggplot(data=auc_df) + 
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

auc_pe = ggarrange(aucplot, mapeplot, nrow=1, ncol = 2, 
                   common.legend = T, legend = "bottom",labels = "AUTO")

ggsave(auc_pe, device = cairo_ps,
       file="report/clinical/images/auc_pe.eps", width = 7, height=5.5)

