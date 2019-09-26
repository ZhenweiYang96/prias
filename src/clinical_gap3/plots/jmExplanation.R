library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")
source("src/clinical_gap3/scheduleCreator.R")

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

pat_data$psa = 2^pat_data$log2psaplus1 - 1

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
  scale_x_continuous(breaks = c(0,1,2,4,6,8,10),
                     limits = c(-0.3, 10)) +
  geom_vline(xintercept = cur_visit_time, linetype="dashed") +
  geom_vline(xintercept = latest_survival_time, color=SUCCESS_COLOR) +
  xlab("Follow-up time (years)")+
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE),
        legend.title = element_blank()) 

blank_common = common + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
        axis.text.y = element_text(size=FONT_SIZE, color=THEME_COLOR))

#c(0,latest_survival_time,cur_visit_time)
A = blank_common + 
  geom_point(aes(x=-50, y=Inf, color='Observed PSA (log scale)')) +
  scale_color_manual(values = THEME_COLOR) + 
  geom_point(aes(x=pat_data$year_visit,y=pat_data$log2psaplus1),
             size=POINT_SIZE, color=THEME_COLOR) +
  geom_line(aes(x=psa_predict_times[psa_predict_times<=cur_visit_time], 
                y=mean_psa[psa_predict_times<=cur_visit_time]), color=THEME_COLOR) + 
  scale_y_continuous(breaks = A_y_breaks, labels = round(A_y_breaks,1)) +
  ylab("PSA (log scale)\nvalue")
#ylab(expression(atop('log'[2]*'(PSA + 1)', 'value')))

B = blank_common + 
  geom_line(aes(x=psa_predict_times[psa_predict_times<=cur_visit_time], 
                y=mean_psa_velocity[psa_predict_times<=cur_visit_time]), 
            color=THEME_COLOR) + 
  scale_y_continuous(breaks = B_y_breaks, labels = round(B_y_breaks,1)) +
  ylab('PSA (log scale)\nInstantaneous\nvelocity')
#ylab(expression(atop('log'[2]*'(PSA + 1)', 'instantaneous velocity')))

C = common + geom_line(aes(x=survival_predict_times, y=mean_cum_risk), color=DANGER_COLOR) +
  geom_ribbon(aes(x=survival_predict_times, ymin=lower_cum_risk,
                  ymax=upper_cum_risk), alpha=0.15, fill=DANGER_COLOR) + 
  scale_y_continuous(breaks = seq(0,1, 0.5), 
                     labels = c("0%", "50%", "100%"), 
                     limits = c(0,1)) + ylab("Cumulative risk\nof reclassification") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color=DANGER_COLOR, size=FONT_SIZE),
        axis.text.y = element_text(color=DANGER_COLOR, size=FONT_SIZE),
        plot.margin = margin(b = 0, unit = "pt"))

D = common +
  geom_label(aes(x=c(0,1,4), y=c(0,0,0), 
                 label = c("Start AS\nGleason\ngrade 1", 
                           "Biopsy\nGleason\ngrade 1",
                           "Current\nVisit")), color='white',
             size= LABEL_SIZE, nudge_x = c(-0.17, 0.17,0),
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
                                heights = c(1.1,1,1.1,0.6), common.legend = T,
                                hjust = -9, vjust = 2,
                                legend = "top", labels = c("A","B", "C", ""))

ggsave(jm_explanation_plot, filename = "report/clinical/images/jmExplanationPlot_113.eps",
       device = cairo_ps, height = 6)


#############
#Now the plot of schedules

E = C + theme(axis.text.x = element_blank())

set.seed(2019)
biopsy_schedules = compareSchedules(object = mvJoint_psa_time_scaled,
                                    patient_data = pat_data,
                                    risk_thresholds = c(0.05, 0.1),
                                    cur_visit_time = cur_visit_time, 
                                    latest_survival_time = latest_survival_time)

schedules = names(biopsy_schedules$schedules)
#check schedule names are in right order
schedules = c("5% Risk", "10% Risk", "Yearly", "Biennially", "PRIAS")

expected_delays = sapply(biopsy_schedules$schedules, "[[", "expected_delay")
total_biopsies = sapply(biopsy_schedules$schedules, "[[", "total_biopsies")

names(expected_delays) = schedules

biopsy_times = do.call('c', lapply(biopsy_schedules$schedules, "[[", "biopsy_times"))
plotDf = data.frame(Schedule=rep(schedules, total_biopsies), 
                    biopsy_times)

plotDf = plotDf[plotDf$Schedule %in% c("5% Risk", "10% Risk", "PRIAS", "Yearly"),]
plotDf$Schedule = droplevels(plotDf$Schedule)

F_ = common +
  geom_line(data = plotDf, aes(x=biopsy_times, 
                               y=as.numeric(plotDf$Schedule),
                               group=as.numeric(plotDf$Schedule)),
            linetype='dotted')+
  geom_label(data = plotDf, 
             aes(x=biopsy_times, y=as.numeric(plotDf$Schedule)),
             label="B",size=POINT_SIZE + 2, fill=DANGER_COLOR, color='white') + 
  ylab("Biopsy schedule") + 
  theme(axis.title.x = element_blank(),
        plot.margin = margin(b = 0, unit = "pt")) +
  scale_y_continuous(breaks = 1:length(levels(plotDf$Schedule)),
                     labels = levels(plotDf$Schedule),
                     limits = c(0.5, length(levels(plotDf$Schedule)) + 0.5)) 

H = ggplot() + geom_bar(aes(x=levels(plotDf$Schedule), 
                            y=12*expected_delays[levels(plotDf$Schedule)]),
                        stat='identity', width=0.5) +  
  ylab("Expected delay (months) in detecting reclassification") + 
  xlab("Biopsy schedule") +
  scale_y_continuous(breaks = seq(0,12, length.out = 7),
                     labels = round(seq(0,12, length.out = 7),1),
                     limits=c(0,12))+
  coord_flip() +
  theme_bw() + 
  theme(text = element_text(size = FONT_SIZE))

demo_pat_plot = ggarrange(E, F_, D, H , ncol=1, nrow=4, align = "v",
                          labels=c("A", "B", "", "C"), heights = c(1,1,0.5,1))

ggsave(demo_pat_plot, filename = "report/clinical/images/demo_pat1.eps",
       device = cairo_ps, height = 6)
