library(ggplot2)

DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
SUCCESS_COLOR = 'forestgreen'
FONT_SIZE = 12
LABEL_SIZE = 2.5

TRUE_TIME_Progression = 4.5
TRUE_TIME_Progression_HT = 0.85

delay_explanation_plot_a = ggplot() + 
  geom_ribbon(aes(x=c(TRUE_TIME_Progression, 5), ymin=-Inf, ymax=Inf), 
                                                  fill=DANGER_COLOR, alpha=0.25) +
  geom_segment(aes(x=0:4, y=rep(-Inf,5), xend=0:4, yend=rep(0.5,5)),
               color=c(WARNING_COLOR, rep(SUCCESS_COLOR, 4)))+
  geom_segment(aes(x=c(TRUE_TIME_Progression,5), y=rep(-Inf,2), xend=c(TRUE_TIME_Progression,5), yend=c(TRUE_TIME_Progression_HT,0.5)),
               color=DANGER_COLOR)+
  geom_label(aes(x=TRUE_TIME_Progression, y=TRUE_TIME_Progression_HT, label = "True time of\nprogression"),
             size=LABEL_SIZE, color=DANGER_COLOR,
             fill='white') +
  geom_label(aes(x=0:5, y=rep(0.5,6),
                 label = c("Start\nsurveillance", 
                           "1st negative\ntest", 
                           "2nd negative\ntest", 
                           "3rd negative\ntest", 
                           "4th negative\ntest",
                           "5th test\nprogression\ndetected")),
             size=LABEL_SIZE, color='white',
             fill=c(WARNING_COLOR, rep(SUCCESS_COLOR,4), DANGER_COLOR)) +
  geom_text(aes(x=0.5 * (TRUE_TIME_Progression + 5), y=0.2, label="6 months delay\n in detecting progression"), size=LABEL_SIZE) + 
  geom_segment(aes(x=TRUE_TIME_Progression, xend = 5, y=0.09, yend = 0.09),
               arrow = arrow(length = unit(1.5,"mm"), ends="both", type="closed"))+
  
  xlab("Follow-up time") + 
  scale_x_continuous(breaks = c(0,1,2,3,TRUE_TIME_Progression,4,5,6), limits=c(-0.2,6.25)) + 
  ylim(0,1) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        title = element_text(size=FONT_SIZE-2),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  ggtitle("    Frequent tests - Shorter delay in detecting progression")


delay_explanation_plot_b = ggplot() + 
  geom_ribbon(aes(x=c(TRUE_TIME_Progression, 6), ymin=-Inf, ymax=Inf), 
                                                  fill=DANGER_COLOR, alpha=0.25) +
  geom_segment(aes(x=c(0,2,4), y=rep(-Inf,3), 
                   xend=c(0,2,4), yend=rep(0.5,3)),
               color=c(WARNING_COLOR, rep(SUCCESS_COLOR, 2)))+
  geom_segment(aes(x=c(TRUE_TIME_Progression,6), y=rep(-Inf,2), 
                   xend=c(TRUE_TIME_Progression,6), yend=c(TRUE_TIME_Progression_HT,0.5)),
               color=DANGER_COLOR)+
  geom_label(aes(x=TRUE_TIME_Progression, y=TRUE_TIME_Progression_HT, label = "True time of\nprogression"),
             size=LABEL_SIZE, color=DANGER_COLOR,
             fill='white') +
  geom_label(aes(x=c(0,2,4,6), y=rep(0.5,4),
                 label = c("Start\nsurveillance", 
                           "1st negative\ntest", 
                           "2nd negative\ntest", 
                           "3rd test\nprogression\ndetected")),
             size=LABEL_SIZE, color='white',
             fill=c(WARNING_COLOR, rep(SUCCESS_COLOR,2), DANGER_COLOR)) +
  geom_text(aes(x=0.5 * (TRUE_TIME_Progression + 6), y=0.2, label="18 months delay\n in detecting progression"), size=LABEL_SIZE) + 
  geom_segment(aes(x=TRUE_TIME_Progression, xend = 6, y=0.09, yend = 0.09),
               arrow = arrow(length = unit(1.5,"mm"), ends="both", type="closed"))+
  xlab("Follow-up time") + 
  scale_x_continuous(breaks = c(0,1,2,3,4, TRUE_TIME_Progression,5,6), 
                     labels = c("Jan 2000","Jan 2001", "Jan 2002",
                                "Jan 2003", "Jan 2004", "Jul 2004", 
                                "Jan 2005", "Jan 2006"), limits=c(-0.2,6.25)) + 
  ylim(0,1) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), 
        title = element_text(size=FONT_SIZE-2),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank())+
  ggtitle("    Infrequent tests - Longer delay in detecting progression")


delay_explanation_plot = ggpubr::ggarrange(delay_explanation_plot_a, delay_explanation_plot_b, 
                                           ncol = 1, nrow=2, align = "v", labels = "AUTO",
                                           heights = c(1, 1.2))

ggsave(delay_explanation_plot, 
       file="report/lastpaper/images/delay_explanation.pdf", 
       device = cairo_pdf, height = 5, width=7/1.333)

ggsave(delay_explanation_plot, 
       file="report/lastpaper/figure1.eps", 
       device = cairo_ps, height = 5, width=7/1.333)
