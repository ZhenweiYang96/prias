load("Rdata/decision_analytic/cleandata.Rdata")
library(JMbayes)
library(splines)

DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
SUCCESS_COLOR = 'forestgreen'

POINT_SIZE = 2
FONT_SIZE = 12
LABEL_SIZE = 2.75

riskGaugeGraph = function(mean_risk_prob, danger_color_threshold = 0.2, gauge_color=DANGER_COLOR){
  
  LABEL_SIZE = 3
  
  risk_label = paste0("\n\n\nReclassification risk\n at current visit: ", round(mean_risk_prob*100), "%")
  
  gauge_ticks_colors = sapply(seq(0,1,0.25), FUN = function(prop){
    if(prop > danger_color_threshold){
      return(DANGER_COLOR)
    }else{
      col=colorRamp(c(SUCCESS_COLOR, WARNING_COLOR, DANGER_COLOR))(prop/danger_color_threshold)
      return(rgb(col[1], col[2], col[3], maxColorValue = 255))
    }
  })
  
  riskGauge = ggplot(data = NULL, 
                     aes(ymax = mean_risk_prob, ymin = 0, xmax = 2, xmin = 1, 
                         fill="Risk")) +
    geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="white", color=gauge_color) +
    geom_rect() +
    geom_segment(aes(x=2.0, xend=2.1, y=0, yend=0), color=gauge_ticks_colors[1])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.25, yend=0.25), color=gauge_ticks_colors[2])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.5, yend=0.5), color=gauge_ticks_colors[3])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.75, yend=0.75), color=gauge_ticks_colors[4])+
    geom_segment(aes(x=2.0, xend=2.1, y=1, yend=1), color=gauge_ticks_colors[5])+
    geom_text(aes(x = 2.55, y = 0, label = "0%"), size=LABEL_SIZE, color=gauge_ticks_colors[1]) +
    geom_text(aes(x = 2.55, y = 0.25, label = "25%"), size=LABEL_SIZE, color=gauge_ticks_colors[2]) +
    geom_text(aes(x = 2.55, y = 0.5, label = "50%"), size=LABEL_SIZE, color=gauge_ticks_colors[3]) +
    geom_text(aes(x = 2.55, y = 0.75, label = "75%"), size=LABEL_SIZE, color=gauge_ticks_colors[4]) +
    geom_text(aes(x = 2.55, y = 1, label = "100%"), size=LABEL_SIZE, color=gauge_ticks_colors[5]) +
    coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2.6)) + ylim(c(0,2)) +
    geom_text(aes(x = 0, y = 0, label = risk_label), 
              color=gauge_color, size=4) +
    geom_segment(aes(x=0, xend=1, y=mean_risk_prob, yend=mean_risk_prob), 
                 color=gauge_color,
                 arrow = arrow(length = unit(0.25,"cm")))+
    geom_point(aes(x=0, y=0), color=gauge_color, size=2)+
    scale_fill_manual("", values=gauge_color)+
    theme_void() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size=FONT_SIZE, color=gauge_color)) +
    guides(fill=FALSE) +
    guides(colour=FALSE)
  return(riskGauge)
}


pat1 = prias_long[prias_long$P_ID==3682,]
pat1 = pat1[!is.na(pat1$psa) & pat1$visitTimeYears<=3,]
pat1$visitTimeYears[nrow(pat1)] = 3
surv1 = 0.95
p1 = ggplot() + 
  geom_point(aes(x=pat1$visitTimeYears, y=pat1$psa), size=POINT_SIZE) + 
  geom_line(aes(x=pat1$visitTimeYears, y=pat1$psa), alpha=0.1) +
  geom_vline(xintercept = 1, color=SUCCESS_COLOR) +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(-0.35,3.25),
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
pat2$visitTimeYears[nrow(pat2)] = 3
surv2 = 0.8
p2 = ggplot() + 
  geom_point(aes(x=pat2$visitTimeYears, y=pat2$psa), size=POINT_SIZE) + 
  geom_line(aes(x=pat2$visitTimeYears, y=pat2$psa), alpha=0.1) +
  geom_vline(xintercept = 1, color=SUCCESS_COLOR) +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE),
        axis.title.x = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "pt"))+
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(-0.35,3.25),
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
  scale_x_continuous(breaks = c(0, 1, 2, 3), limits = c(-0.35,3.25),
                     labels = c("0", "1","2", "3"))


psa_plot = ggpubr::ggarrange(p1, p2, p3, 
                             align = "v", labels = c("A", "B", ""),
                             ncol=1, nrow=3, heights = c(1.1, 1, 0.25))

risk_plot = ggpubr::ggarrange(ggplot() + theme_void(),
                              riskGaugeGraph(mean_risk_prob = 0.05, gauge_color = SUCCESS_COLOR),
                              riskGaugeGraph(mean_risk_prob = 0.2, gauge_color = DANGER_COLOR),
                              ggplot() + theme_void(),
                              ncol = 1, nrow = 4, align = "v", heights = c(0.35,1,1, 0.15))

final_plot = ggpubr::ggarrange(psa_plot, risk_plot, widths = c(2,1))

print(final_plot)

ggsave(final_plot, device = cairo_ps,
       height = 5.5,
       filename = "report/clinical/images/riskBasedExample.eps")
