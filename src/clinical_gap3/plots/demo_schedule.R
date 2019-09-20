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

accuracy = 50
survival_predict_times = seq(latest_survival_time, MAX_FOLLOW_UP, length.out = accuracy)

exp_fut = getExpectedFutureOutcomes(mvJoint_psa_time_scaled, pat_data, latest_survival_time, Inf,
                                    survival_predict_times, psaDist = "Tdist")

mean_cum_risk = 1-c(1, rowMeans(exp_fut$predicted_surv_prob))
lower_cum_risk = 1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.025))
upper_cum_risk = 1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.975))

common = ggplot() +  
  scale_x_continuous(breaks = xbreaks, labels=xlabs,
                     limits = c(0, max_follow_up), 
                     minor_breaks = seq(0, max_follow_up, 1)) +
  xlab("Follow-up time (years)")+
  geom_vline(xintercept = latest_survival_time) + 
  geom_vline(xintercept = max_psa_time, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE)) 

blank_common = common + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank())

B = blank_common + geom_line(aes(x=survival_predict_times, y=mean_cum_risk), color=DANGER_COLOR) +
  geom_ribbon(aes(x=survival_predict_times, ymin=lower_cum_risk,
                  ymax=upper_cum_risk), alpha=0.15, fill=DANGER_COLOR) + 
  scale_y_continuous(breaks = seq(0,1, 0.5), 
                     labels = c("0%", "50%", "100%"), 
                     limits = c(0,1)) + ylab("Cumulative risk\nof reclassification")+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color=DANGER_COLOR, size=FONT_SIZE),
        axis.text.y = element_text(color=DANGER_COLOR, size=FONT_SIZE))

biopsy_schedules = compareSchedules(mvJoint_psa_time_scaled, pat_data,cur_visit_time = max_psa_time, 
                                    latest_survival_time = latest_survival_time)

schedules = c("5% Risk", "10% Risk", "15% Risk", 
              "PRIAS", "Yearly", "Every 2 Years")
expected_delays = sapply(biopsy_schedules$schedules, "[[", "expected_delay")
total_biopsies = sapply(biopsy_schedules$schedules, "[[", "total_biopsies")

biopsy_times = do.call('c', lapply(biopsy_schedules$schedules, "[[", "biopsy_times"))
plotDf = data.frame(Schedule=rep(schedules, total_biopsies), 
                    biopsy_times)

plotDf = plotDf[plotDf$Schedule %in% c("5% Risk", "10% Risk", "PRIAS", "Yearly"),]
plotDf$Schedule = droplevels(plotDf$Schedule)

C = common +
  geom_line(data = plotDf, aes(x=biopsy_times, 
                               y=as.numeric(plotDf$Schedule),
                               group=as.numeric(plotDf$Schedule)),
            linetype='dotted')+
  geom_label(data = plotDf, 
             aes(x=biopsy_times, y=as.numeric(plotDf$Schedule)),
             label="B",size=POINT_SIZE + 2, fill='black', color='white') + 
  ylab("Biopsy schedule") + 
  scale_y_continuous(breaks = 1:length(levels(plotDf$Schedule)),
                     labels = levels(plotDf$Schedule),
                     limits = c(0.5, length(levels(plotDf$Schedule)) + 0.5)) 


D = ggplot() + geom_bar(aes(x=levels(plotDf$Schedule), 
                            y=12*expected_delays[c("Risk: 0.1", "Risk: 0.05", "prias", "annual")]),
                        stat='identity', width=0.5) +  
  ylab("Expected delay (months)\n in detecting Gleason \u2265 7") + 
  xlab("Biopsy schedule") +
  scale_y_continuous(breaks = seq(0,12, length.out = 4),
                     labels = round(seq(0,12, length.out = 4),1),
                     limits=c(0,12))+
  coord_flip() +
  theme_bw() + 
  theme(text = element_text(size = FONT_SIZE))

E = ggplot() + geom_bar(aes(x=levels(plotDf$Schedule), 
                            y=total_biopsies[c("Risk: 0.1", "Risk: 0.05", "prias", "annual")]),
                        stat='identity', width=0.5) +  
  ylab("Total biopsies scheduled") + 
  xlab("Biopsy schedule") + 
  theme_bw() + 
  theme(text = element_text(size = FONT_SIZE),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()

upper_plot = ggarrange(B, C, ncol=1, nrow=2, align = "v",
                       labels="AUTO", heights = c(1,1.35))
lower_plot = ggarrange(D, E, ncol=2, nrow=1, 
                       align="h", widths = c(1.25,1))

demo_pat_plot = ggarrange(upper_plot, lower_plot, ncol=1, nrow=2, heights = c(3,1))


ggsave(demo_pat_plot, filename = "report/clinical/images/demo_pat1.eps",
       device = cairo_ps, height = 9, width = 7)
