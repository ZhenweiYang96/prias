library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")

SUCCESS_COLOR = 'forestgreen'
DANGER_COLOR = 'red'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10

DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
SUCCESS_COLOR = 'forestgreen'

POINT_SIZE = 2
FONT_SIZE = 11
LABEL_SIZE = 2.75

pat_data = prias_long_final[prias_long_final$P_ID==113,]
set.seed(2019)
pat_data$log2psaplus1 = pat_data$log2psaplus1 - 2

pat_data$log2psaplus1[pat_data$year_visit > 2.5] = pat_data$log2psaplus1[pat_data$year_visit > 2.5] + 
  rnorm(n = length(pat_data$log2psaplus1[pat_data$year_visit > 2.5]), 2.25, 1)

latest_survival_time = 1

cur_visit_time = max(pat_data$year_visit)

accuracy = 100
psa_predict_times = seq(0, MAX_FOLLOW_UP, length.out = accuracy)
survival_predict_times = seq(latest_survival_time, MAX_FOLLOW_UP, length.out = accuracy)

exp_fut = getExpectedFutureOutcomes(mvJoint_psa_time_scaled, pat_data, latest_survival_time, Inf,
                                    survival_predict_times, psa_predict_times, psaDist = "Tdist")
mean_psa = rowMeans(exp_fut$predicted_psa)
mean_psa_velocity = rowMeans(exp_fut$predicted_psa_slope)
mean_cum_risk = 1-c(1, rowMeans(exp_fut$predicted_surv_prob))
lower_cum_risk = 1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.025))
upper_cum_risk = 1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.975))

A_y_breaks = seq(min(pat_data$log2psaplus1, na.rm = T), 
                 max(pat_data$log2psaplus1, na.rm = T), 
                 length.out = 3)

B_y_breaks = seq(min(mean_psa_velocity, na.rm = T), 
                 max(mean_psa_velocity, na.rm = T), 
                 length.out = 3)

common = ggplot() +  
  scale_x_continuous(breaks = 0:4, labels=0:4,
                     limits = c(-0.25, 4+0.25), 
                     minor_breaks = seq(0, 4, 1)) +
  xlab("Follow-up time (years)")+
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE)) 

blank_common = common + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
        axis.text.y = element_text(size=FONT_SIZE, color=THEME_COLOR))

#c(0,latest_survival_time,cur_visit_time)
A = blank_common + 
  geom_point(aes(x=pat_data$year_visit,y=pat_data$log2psaplus1),
             size=POINT_SIZE, color=THEME_COLOR) +
  geom_line(aes(x=psa_predict_times, y=mean_psa), color=THEME_COLOR) + 
  scale_y_continuous(breaks = A_y_breaks, labels = round(A_y_breaks,1)) +
  ylab("PSA (log scale)\nvalue")
#ylab(expression(atop('log'[2]*'(PSA + 1)', 'value')))

B = blank_common + 
  geom_line(aes(x=psa_predict_times, y=mean_psa_velocity), color=THEME_COLOR) + 
  scale_y_continuous(breaks = B_y_breaks, labels = round(B_y_breaks,1)) +
  ylab('PSA (log scale)\nInstantaneous\nvelocity')
#ylab(expression(atop('log'[2]*'(PSA + 1)', 'instantaneous velocity')))

C = common + geom_line(aes(x=survival_predict_times, y=mean_cum_risk), color=DANGER_COLOR) +
  geom_ribbon(aes(x=survival_predict_times, ymin=lower_cum_risk,
                  ymax=upper_cum_risk), alpha=0.15, fill=DANGER_COLOR) + 
  geom_vline(xintercept = latest_survival_time, color=SUCCESS_COLOR) +
  geom_point(aes(x=-50, y=Inf, color='Observed PSA (log scale)')) +
  scale_color_manual(values = THEME_COLOR) + 
  scale_y_continuous(breaks = seq(0,1, 0.5), 
                     labels = c("0%", "50%", "100%"), 
                     limits = c(0,1)) + ylab("Cumulative risk\nof reclassification") +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color=DANGER_COLOR, size=FONT_SIZE),
        axis.text.y = element_text(color=DANGER_COLOR, size=FONT_SIZE),
        plot.margin = margin(b = 0, unit = "pt"))

D = common +
  geom_label(aes(x=c(0,1,4), y=c(0,0,0), 
                 label = c("Start AS\nGleason\ngrade 1", 
                           "Biopsy\nGleason\ngrade 1",
                           "Current\nVisit")), color='white',
             size= LABEL_SIZE,
             fill=c(WARNING_COLOR, SUCCESS_COLOR, 'black')) +
  theme(text = element_text(size = FONT_SIZE),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) + 
  ylim(-0.25,0.25)

jm_explanation_plot = ggarrange(A,B, C, D, ncol=1, nrow=4, align = "v",
                                heights = c(1.1,1,1.1,0.55), common.legend = T,
                                hjust = -9, vjust = 2,
                                legend = "top", labels = c("A","B", "C", ""))

ggsave(jm_explanation_plot, filename = "report/clinical/images/jmExplanationPlot_113.eps",
       device = cairo_ps, height = 6)



