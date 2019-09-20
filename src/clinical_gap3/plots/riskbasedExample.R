load("Rdata/decision_analytic/cleandata.Rdata")
library(JMbayes)
library(splines)

DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
SUCCESS_COLOR = 'forestgreen'

POINT_SIZE = 2
FONT_SIZE = 12
LABEL_SIZE = 2.75

pat1 = prias_long[prias_long$P_ID==3682,]
pat1 = pat1[!is.na(pat1$psa) & pat1$visitTimeYears<=3,]
surv1 = 0.95
p1 = ggplot() + geom_col(aes(x=c(3,3), y=15 * c(1-surv1, surv1)), 
                         fill=c(SUCCESS_COLOR, "white"), 
                         color=SUCCESS_COLOR, width = 0.4) + 
  geom_text(aes(x=3, y=2.5, label='Reclassification\nRisk 5%'),
            color=SUCCESS_COLOR, size=LABEL_SIZE) +
  geom_point(aes(x=pat1$visitTimeYears, y=pat1$psa), size=POINT_SIZE) + 
  geom_line(aes(x=pat1$visitTimeYears, y=pat1$psa), alpha=0.1) +
  geom_vline(xintercept = 1, color=SUCCESS_COLOR) +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(-0.25,3.5),
                     labels = c("0", "1","2", "3")) + 
  ylim(0,15) +
  ylab("PSA (ng/mL)") + xlab("Follow-up time (years)") +
  ggtitle("Should a biopsy be conducted at current visit?")

set.seed(100)
pat2 = prias_long[prias_long$P_ID==1931,]
pat2 = pat2[!is.na(pat2$psa) & pat2$visitTimeYears<=3,]
pat2$psa[pat2$visitTimeYears>1] = pat2$psa[pat2$visitTimeYears>1] + 
  2*pat2$visitTimeYears[pat2$visitTimeYears>1] + 
  rnorm(n=sum(pat2$visitTimeYears>2), mean = 0, sd=0.18)
surv2 = 0.8
p2 = ggplot() + geom_col(aes(x=c(3,3), y=15 * c(1-surv2, surv2)), 
                         fill=c(DANGER_COLOR, "white"),  color=DANGER_COLOR, width = 0.4) + 
  geom_text(aes(x=3, y=4.5, label='Reclassification\nRisk 20%'), color=DANGER_COLOR, size=LABEL_SIZE) +
  geom_point(aes(x=pat2$visitTimeYears, y=pat2$psa), size=POINT_SIZE) + 
  geom_line(aes(x=pat2$visitTimeYears, y=pat2$psa), alpha=0.1) +
  geom_vline(xintercept = 1, color=SUCCESS_COLOR) +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE),
        axis.title.x = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "pt"))+
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(-0.25,3.5),
                     labels = c("0", "1","2", "3")) + 
  ylim(0,15) +
  ylab("PSA (ng/mL)") + xlab("Follow-up time (years)")

p3 = ggplot() + 
  geom_label(aes(x=c(0,1,3), y=c(0,0,0), 
                 label = c("Start AS\nGleason grade 1", 
                           "Biopsy\nGleason grade 1",
                           "Current\nVisit")), color='white',
             size= LABEL_SIZE,
             fill=c(WARNING_COLOR, SUCCESS_COLOR, 'black')) +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "pt")) + 
  xlab("Follow-up time (years)") + ylim(-0.25,0.25) + 
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(-0.25,3.5),
                     labels = c("0", "1","2", "3"))


final_plot = ggpubr::ggarrange(p1, p2, p3, 
                  align = "v", labels = c("A", "B", ""),
                  ncol=1, nrow=3, heights = c(1.1, 1, 0.25))
print(final_plot)

ggsave(final_plot, device = cairo_ps,
       height = 5.5,
       filename = "report/clinical/images/riskBasedExample.eps")
