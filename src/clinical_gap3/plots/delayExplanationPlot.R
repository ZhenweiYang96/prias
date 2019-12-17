library(ggplot2)

DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
SUCCESS_COLOR = 'forestgreen'
FONT_SIZE = 12
LABEL_SIZE = 2.75

TRUE_TIME_GS7 = 3.5
TRUE_TIME_GS7_HT = 0.85

delay_explanation_plot_a = ggplot() + geom_ribbon(aes(x=c(TRUE_TIME_GS7, 4), ymin=-Inf, ymax=Inf), 
                                                  fill=DANGER_COLOR, alpha=0.25) +
  geom_segment(aes(x=0:3, y=rep(-Inf,4), xend=0:3, yend=rep(0.5,4)),
               color=c(WARNING_COLOR, rep(SUCCESS_COLOR, 3)))+
  geom_segment(aes(x=c(TRUE_TIME_GS7,4), y=rep(-Inf,2), xend=c(TRUE_TIME_GS7,4), yend=c(TRUE_TIME_GS7_HT,0.5)),
               color=DANGER_COLOR)+
  geom_label(aes(x=TRUE_TIME_GS7, y=TRUE_TIME_GS7_HT, label = "True time of\nGleason grade \u2265 2"),
             size=LABEL_SIZE, color=DANGER_COLOR,
             fill='white') +
  geom_label(aes(x=0:4, y=rep(0.5,5),
                 label = c("Start AS\nGleason\ngrade 1", 
                           "1st Biopsy\nGleason\ngrade 1", 
                           "2nd Biopsy\nGleason\ngrade 1", 
                           "3rd Biopsy\nGleason\ngrade 1", 
                           "4th Biopsy\nGleason\ngrade \u2265 2")),
             size=LABEL_SIZE, color='white',
             fill=c(WARNING_COLOR, rep(SUCCESS_COLOR,3), DANGER_COLOR)) +
  geom_text(aes(x=0.5 * (TRUE_TIME_GS7 + 4), y=0.2, label="6 months delay\n in detecting upgrading"), size=LABEL_SIZE) + 
  geom_segment(aes(x=TRUE_TIME_GS7, xend = 4, y=0.09, yend = 0.09),
               arrow = arrow(length = unit(2.5,"mm"), ends="both", type="closed"))+
  
  xlab("Time of biopsy visits") + 
  scale_x_continuous(breaks = c(0,1,2,3,TRUE_TIME_GS7,4), 
                     labels = c("Jan 2005","Jan 2006", "Jan 2007", "Jan 2008",
                                "Jul 2008", "Jan 2009"), limits=c(-0.2,5.25)) + 
  ylim(0,1) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  ggtitle("    Biopsy every year")


delay_explanation_plot_b = ggplot() + geom_ribbon(aes(x=c(TRUE_TIME_GS7, 5), ymin=-Inf, ymax=Inf), 
                                                  fill=DANGER_COLOR, alpha=0.25) +
  geom_segment(aes(x=c(0,1,3), y=rep(-Inf,3), xend=c(0,1,3), yend=rep(0.5,3)),
               color=c(WARNING_COLOR, rep(SUCCESS_COLOR, 2)))+
  geom_segment(aes(x=c(TRUE_TIME_GS7,5), y=rep(-Inf,2), xend=c(TRUE_TIME_GS7,5), yend=c(TRUE_TIME_GS7_HT,0.5)),
               color=DANGER_COLOR)+
  geom_label(aes(x=TRUE_TIME_GS7, y=TRUE_TIME_GS7_HT, label = "True time of\nGleason grade \u2265 2"),
             size=LABEL_SIZE, color=DANGER_COLOR,
             fill='white') +
  geom_label(aes(x=c(0,1,3,5), y=rep(0.5,4),
                 label = c("Start AS\nGleason\ngrade 1", 
                           "1st Biopsy\nGleason\ngrade 1", 
                           "2nd Biopsy\nGleason\ngrade 1", 
                           "3rd Biopsy\nGleason\ngrade \u2265 2")),
             size=LABEL_SIZE, color='white',
             fill=c(WARNING_COLOR, rep(SUCCESS_COLOR,2), DANGER_COLOR)) +
  geom_text(aes(x=0.5 * (TRUE_TIME_GS7 + 5), y=0.2, label="18 months delay\n in detecting upgrading"), size=LABEL_SIZE) + 
  geom_segment(aes(x=TRUE_TIME_GS7, xend = 5, y=0.09, yend = 0.09),
               arrow = arrow(length = unit(2.5,"mm"), ends="both", type="closed"))+
  xlab("Time of biopsy visits") + 
  scale_x_continuous(breaks = c(0,1,2,3,TRUE_TIME_GS7,4, 5), 
                     labels = c("Jan 2005","Jan 2006", "Jan 2007","Jan 2008",
                                "Jul 2008", "Jan 2009", "Jan 2010"), limits=c(-0.2,5.25)) + 
  ylim(0,1) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank())+
  ggtitle("    Biopsy every 2 years")


delay_explanation_plot = ggpubr::ggarrange(delay_explanation_plot_a, delay_explanation_plot_b, 
                                           ncol = 1, nrow=2, align = "v", labels = "AUTO",
                                           heights = c(1, 1.2))

ggsave(delay_explanation_plot, 
       file="report/clinical/images/delay_explanation.eps", 
       device = cairo_ps, height = 5.5, width = 6)
