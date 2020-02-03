library(ggplot2)
library(ggpubr)
FONT_SIZE = 12

seeds = c(2001:2100)

biopsyDf_summary = do.call('rbind', lapply(seeds, FUN = function(seed){
  load(paste0("Rdata/lastpaper/simulation/combined_results_non_auto/seed_", seed, ".Rdata"))
  
  return(biopsyDf_summary)
}))

biopsyDf_summary = droplevels(biopsyDf_summary[biopsyDf_summary$schedule %in%
                                      c("Annual", "PRIAS", "Risk: 10%", "Risk: 5%", "Risk: Auto (0.5)"),])

levels(biopsyDf_summary$schedule)[5] = "Risk:\u03BA*(v)"
# 
# biopsyDf_summary$nb[biopsyDf_summary$delay<0] = biopsyDf_summary$nb[biopsyDf_summary$delay<0] + 1
# biopsyDf_summary$delay[biopsyDf_summary$delay<0] = 10 - biopsyDf_summary$progression_time[biopsyDf_summary$delay<0]

by(INDICES = biopsyDf_summary$schedule[biopsyDf_summary$progressed==0],
   data= biopsyDf_summary$nb[biopsyDf_summary$progressed==0], summary)

by(INDICES = biopsyDf_summary$schedule[biopsyDf_summary$progressed==1],
   data= biopsyDf_summary$nb[biopsyDf_summary$progressed==1], summary)

by(INDICES = biopsyDf_summary$schedule[biopsyDf_summary$progressed==1],
   data= biopsyDf_summary$delay[biopsyDf_summary$progressed==1], summary)


a = ggplot(biopsyDf_summary[biopsyDf_summary$progressed==1,]) +
  geom_boxplot(aes(x=schedule, y=nb), outlier.shape = NA) +
  theme_bw() + 
  scale_y_continuous(breaks = c(1:10), limits = c(1,10)) +
  theme(text=element_text(size=FONT_SIZE),
        title = element_text(size=FONT_SIZE-1)) +
  xlab("Schedule") + ylab("Number of biopsies") +
  coord_flip() +
  ggtitle("Progressing patients (50%)")

b = ggplot(biopsyDf_summary[biopsyDf_summary$progressed==1,]) +
  geom_boxplot(aes(x=schedule, y=delay), outlier.shape = NA)+
  theme_bw() + 
  scale_y_continuous(breaks = 0:4, limits = c(0,4)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        text=element_text(size=FONT_SIZE)) +
  xlab("Schedule") + ylab("Time delay in detection\nof progression (years)")+
  coord_flip()

c = ggplot(biopsyDf_summary[biopsyDf_summary$progressed==0,]) +
  geom_boxplot(aes(x=schedule, y=nb), outlier.shape = NA) +
  theme_bw() +
  scale_y_continuous(breaks = c(1:10), limits = c(1,10)) +
  theme(text=element_text(size=FONT_SIZE),
        title = element_text(size=FONT_SIZE-1)) +
  xlab("Schedule") + ylab("Number of biopsies") +
  coord_flip() +
  ggtitle("Non-progressing patients (50%)")

d = ggplot() +
  geom_boxplot(data=biopsyDf_summary[biopsyDf_summary$progressed==0,],
               aes(x=schedule, y=-1), outlier.shape = NA)+
  geom_text(aes(x="Risk: 10%", y=2), label="Time delay not available for\nnon-progressing patients")+
  theme_bw() + 
  scale_y_continuous(breaks = 0:4, limits = c(0,4)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        text=element_text(size=FONT_SIZE)) +
  xlab("Schedule") + ylab("Time delay in detection\nof progression (years)")+ 
  coord_flip()

upper_plot = ggarrange(a,b, nrow=1, ncol=2, align = "h", widths = c(1.225,1))
lower_plot = ggarrange(c,d, nrow=1, ncol=2, align = "h", widths = c(1.225,1))

final_plot = ggpubr::ggarrange(upper_plot, lower_plot, 
                  nrow=2, ncol=1, align = "v", labels = "AUTO")

print(final_plot)
ggsave(final_plot, filename = "report/lastpaper/images/simulation_boxplot.eps",
       device = cairo_ps,  height = 6.5, width=6.5)

