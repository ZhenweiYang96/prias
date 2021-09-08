library(ggplot2)
library(ggpubr)

load(file = "Rdata/lastpaper/simulation/sim_res.Rdata")

sim_res = droplevels(sim_res[sim_res$schedule %in%
                                                 c("Annual", "PRIAS", "Risk: 10%", "Risk: Auto (Inf)", "Risk: Auto (0.75)"),])
levels(sim_res$schedule)

levels(sim_res$schedule)[3] = "\u03BA=10%"
levels(sim_res$schedule)[4] = "\u03BA*{v | E(D) \u2264 0.75}"
levels(sim_res$schedule)[5] = "\u03BA*(v)"
# 
# sim_res$nb[sim_res$delay<0] = sim_res$nb[sim_res$delay<0] + 1
# sim_res$delay[sim_res$delay<0] = 10 - sim_res$progression_time[sim_res$delay<0]

by(INDICES = sim_res$schedule[sim_res$progressed==0],
   data= sim_res$nb[sim_res$progressed==0], summary)

by(INDICES = sim_res$schedule[sim_res$progressed==1],
   data= sim_res$nb[sim_res$progressed==1], summary)

by(INDICES = sim_res$schedule[sim_res$progressed==1],
   data= sim_res$delay[sim_res$progressed==1], summary)


FONT_SIZE = 12
POINT_SIZE = 1.5
THEME_COLOR = 'dodgerblue4'
WARNING_COLOR = 'darkorange'

a = ggplot(sim_res[sim_res$progressed==1,], aes(x=schedule, y=nb)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom="point", size=POINT_SIZE, color=WARNING_COLOR) +
  theme_bw() + 
  scale_y_continuous(breaks = c(1,4,7,10), 
                     limits = c(1,10), minor_breaks = 1:10) +
  theme(text=element_text(size=FONT_SIZE),
        title = element_text(size=FONT_SIZE-2),
        axis.title = element_text(size=FONT_SIZE)) +
  xlab("Schedule") + ylab("Number of biopsies") +
  coord_flip() +
  ggtitle("Progressing patients (50%)")

b = ggplot(sim_res[sim_res$progressed==1,], aes(x=schedule, y=delay)) +
  geom_boxplot(outlier.shape = NA)+
  stat_summary(fun.y=mean, geom="point", size=POINT_SIZE, color=WARNING_COLOR) +
  theme_bw() + 
  scale_y_continuous(breaks = 0:3, limits = c(0,3)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        text=element_text(size=FONT_SIZE)) +
  xlab("Schedule") + ylab("Time delay in detecting\nprogression (years)")+
  coord_flip()

c = ggplot(sim_res[sim_res$progressed==0,], aes(x=schedule, y=nb)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom="point", size=POINT_SIZE, color=WARNING_COLOR) +
  theme_bw() +
  scale_y_continuous(breaks = c(1,4,7,10), limits = c(1,10),
                     minor_breaks = 1:10) +
  theme(text=element_text(size=FONT_SIZE),
        title = element_text(size=FONT_SIZE-2),
        axis.title = element_text(size=FONT_SIZE)) +
  xlab("Schedule") + ylab("Number of biopsies") +
  coord_flip() +
  ggtitle("Non-progressing patients (50%)")

d = ggplot() + geom_text(aes(x=1, y=1.5), label="Time delay not available for\nnon-progressing patients")+
  theme_bw() + 
  scale_y_continuous(breaks = 0:3, limits = c(0,3)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        text=element_text(size=FONT_SIZE), axis.ticks.y=element_blank()) +
  xlab("Schedule") + ylab("Time delay in detecting\nprogression (years)")+ 
  coord_flip()

upper_plot = ggarrange(a,b, nrow=1, ncol=2, align = "h", widths = c(1.5,1))
lower_plot = ggarrange(c,d, nrow=1, ncol=2, align = "h", widths = c(1.5,1))

final_plot = ggpubr::ggarrange(upper_plot, lower_plot, 
                  nrow=2, ncol=1, align = "v", labels = "AUTO")

print(final_plot)
ggsave(final_plot, filename = "report/lastpaper/images/simulation_boxplot.pdf",
       device = cairo_pdf,  height = 6.5, width=6.5)
ggsave(final_plot, filename = "report/lastpaper/figure6.eps",
       device = cairo_ps,  height = 6.5, width=6.5)

