load("Rdata/decision_analytic/cleandata.Rdata")
source("src/mdp/common/prediction.R")
load("Rdata/decision_analytic/PSA_Only/mvJoint_psa.Rdata")
library(JMbayes)
library(splines)

POINT_SIZE = 3
FONT_SIZE = 18

pat1 = prias_long[prias_long$P_ID==3682,]
pat1 = pat1[!is.na(pat1$psa) & pat1$visitTimeYears<=3,]
surv1 = survfitJM(mvJoint_psa, newdata = pat1, idVar = 'P_ID', 
          survTimes = 3, last.time = 1)[[1]][[1]][,'Mean']
surv1 = 0.95
p1 = ggplot() + geom_col(aes(x=c(3,3), y=15 * c(1-surv1, surv1)), 
                    fill=c("forestgreen", "white"), 
                    color="forestgreen", width = 0.5) + 
  geom_text(aes(x=3, y=2.5, label='Risk = 5%'),
            color='forestgreen', size=4) +
  geom_point(aes(x=pat1$visitTimeYears, y=pat1$psa), size=POINT_SIZE) + 
  geom_line(aes(x=pat1$visitTimeYears, y=pat1$psa), alpha=0.1) +
  geom_vline(xintercept = 1) +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(0,3)) + 
  ylim(0,15) +
  ylab("PSA (ng/mL)") + xlab("Follow-up time (years)")

set.seed(100)
pat2 = prias_long[prias_long$P_ID==1931,]
pat2 = pat2[!is.na(pat2$psa) & pat2$visitTimeYears<=3,]
pat2$psa[pat2$visitTimeYears>1] = pat2$psa[pat2$visitTimeYears>1] + 
  2*pat2$visitTimeYears[pat2$visitTimeYears>1] + 
  rnorm(n=sum(pat2$visitTimeYears>2), mean = 0, sd=0.18)
surv2 = survfitJM(mvJoint_psa, newdata = pat2, idVar = 'P_ID', 
                  survTimes = 3, last.time = 1)[[1]][[1]][,'Mean']
surv2 = 0.8
p2 = ggplot() + geom_col(aes(x=c(3,3), y=15 * c(1-surv2, surv2)), 
                    fill=c("red", "white"),  color="red", width = 0.5) + 
  geom_text(aes(x=3, y=3.5, label='Risk = 20%'), color='red', size=4) +
  geom_point(aes(x=pat2$visitTimeYears, y=pat2$psa), size=POINT_SIZE) + 
  geom_line(aes(x=pat2$visitTimeYears, y=pat2$psa), alpha=0.1) +
  geom_vline(xintercept = 1) +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE))+
  scale_x_continuous(breaks = c(0,1,3),
                     labels = c("t=0\n(AS started)",
                                "t=1\n(Latest negative biopsy)",
                                "t=3\n(Current visit)")) + 
  ylim(0,15) +
  ylab("PSA (ng/mL)") + xlab("Follow-up time (years)")

final_plot = ggpubr::ggarrange(p1, p2, labels = "AUTO",
                  ncol=1, nrow=2,heights = c(1, 1.2))

ggsave(final_plot, device = cairo_ps,
       height = 8,
       filename = "report/clinical/images/riskBasedExample.eps")