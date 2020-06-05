library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")
source("src/clinical_gap3/fixedAndPRIASSchedule.R")
source("src/lastpaper/pers_schedule_api.R")

SUCCESS_COLOR = 'forestgreen'
DANGER_COLOR = 'red'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 6

DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
SUCCESS_COLOR = 'forestgreen'

POINT_SIZE = 2
FONT_SIZE = 11
LABEL_SIZE = 2.75

pat_data = prias_long_final[prias_long_final$P_ID==3417 & prias_long_final$year_visit<=2,]
set.seed(2019)

latest_survival_time = 1
pat_data$year_visit[nrow(pat_data)] = 2
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
  scale_x_continuous(breaks = 0:MAX_FOLLOW_UP,
                     limits = c(-0.3, MAX_FOLLOW_UP)) +
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
  geom_point(aes(x=-50, y=Inf, color='Observed PSA (log scale of ng/mL)')) +
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
                     limits = c(0,1)) + ylab("Cause-specific\ncumulative upgrading-risk (%)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color=DANGER_COLOR, size=FONT_SIZE),
        axis.text.y = element_text(color=DANGER_COLOR, size=FONT_SIZE),
        plot.margin = margin(b = 0, unit = "pt"))

D = common +
  geom_label(aes(x=c(0,1,cur_visit_time), y=c(0,0,0), 
                 label = c("Start AS\nGleason\ngrade 1", 
                           "Biopsy\nGleason\ngrade 1",
                           "Current\nVisit")), color='white',
             size= LABEL_SIZE, nudge_x = c(-0.1, 0, 0),
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
print(jm_explanation_plot)

ggsave(jm_explanation_plot, filename = "report/clinical/images/jmExplanationPlot_113.eps",
       device = cairo_ps, height = 7, width=6)


################
#Now the schedules

E = C + theme(axis.text.x = element_blank())

set.seed(2019)
M=750

planned_prias_schedule = getPRIASSchedule(object = mvJoint_psa_time_scaled,
                                          patient_data = pat_data, 
                                          cur_visit_time = cur_visit_time, 
                                          latest_survival_time = latest_survival_time,
                                          M=M, horizon = MAX_FOLLOW_UP)
if(max(planned_prias_schedule)<MAX_FOLLOW_UP){
  planned_prias_schedule = c(planned_prias_schedule, MAX_FOLLOW_UP)
}

prias_schedule = testScheduleConsequences.mvJMbayes(object = mvJoint_psa_time_scaled,
                                                    newdata = pat_data, idVar = "P_ID", 
                                                    last_test_time = latest_survival_time,
                                                    planned_test_schedule = planned_prias_schedule,
                                                    seed = 2019, M = M)

all_risk_schedules = personalizedSchedule.mvJMbayes(object = mvJoint_psa_time_scaled,
                                                    newdata = pat_data, idVar = "P_ID", 
                                                    last_test_time = latest_survival_time,
                                                    fixed_grid_visits = seq(cur_visit_time, MAX_FOLLOW_UP, 0.5),
                                                    gap = 1, seed = 2019, M = M)

risk_5_perc_schedule = all_risk_schedules$all_schedules$`Risk: 0.05`
risk_10_perc_schedule = all_risk_schedules$all_schedules$`Risk: 0.11`
annual_schedule = all_risk_schedules$all_schedules$`Risk: 0`

schedule_df = data.frame(Schedule=character(), number=numeric(),
                         biopsy_times=vector(),
                         expected_num_tests=numeric(),
                         expected_detection_delay=numeric())
schedule_df = rbind(schedule_df, data.frame(Schedule="Yearly", 
                                            number=1,
                                            biopsy_times=annual_schedule$planned_test_schedule,
                                            expected_num_tests=annual_schedule$expected_num_tests,
                                            expected_detection_delay=annual_schedule$expected_detection_delay))
schedule_df = rbind(schedule_df, data.frame(Schedule="PRIAS", 
                                            number=2,
                                            biopsy_times=prias_schedule$planned_test_schedule,
                                            expected_num_tests=prias_schedule$expected_num_tests,
                                            expected_detection_delay=prias_schedule$expected_detection_delay))
schedule_df = rbind(schedule_df, data.frame(Schedule="10% Risk", 
                                            number=3,
                                            biopsy_times=risk_10_perc_schedule$planned_test_schedule,
                                            expected_num_tests=risk_10_perc_schedule$expected_num_tests,
                                            expected_detection_delay=risk_10_perc_schedule$expected_detection_delay))
schedule_df = rbind(schedule_df, data.frame(Schedule="5% Risk", 
                                            number=4,
                                            biopsy_times=risk_5_perc_schedule$planned_test_schedule,
                                            expected_num_tests=risk_5_perc_schedule$expected_num_tests,
                                            expected_detection_delay=risk_5_perc_schedule$expected_detection_delay))

F_ = common +
  geom_line(data = schedule_df, aes(x=biopsy_times, 
                               y=as.numeric(schedule_df$Schedule),
                               group=as.numeric(schedule_df$Schedule)),
            linetype='dotted')+
  geom_label(data = schedule_df, 
             aes(x=biopsy_times, y=as.numeric(schedule_df$Schedule)),
             label="B",size=POINT_SIZE + 2, fill=DANGER_COLOR, color='white') + 
  ylab("Biopsy schedule") + 
  theme(axis.title.x = element_blank(),
        plot.margin = margin(b = 0, unit = "pt")) +
  scale_y_continuous(breaks = 1:length(levels(schedule_df$Schedule)),
                     labels = levels(schedule_df$Schedule),
                     limits = c(0.5, length(levels(schedule_df$Schedule)) + 0.5)) 


consequences_df = schedule_df[!duplicated(schedule_df$Schedule),]
max_delay_limit = 3
H = ggplot() + geom_col(aes(x=rep(consequences_df$Schedule,2), 
                                     y=c(consequences_df$expected_detection_delay, 
                                         max_delay_limit - consequences_df$expected_detection_delay)),
                                 color='black', fill=c(rep(c('darkgrey','white'), nrow(consequences_df))),
                                 width=0.5)+
  ylab("Expected time delay (years)\nin detecting upgrading") + 
  xlab("Biopsy Schedule") +
  scale_y_continuous(breaks = seq(0, max_delay_limit, by = 0.5),
                     labels = seq(0, max_delay_limit,  by = 0.5),
                     limits= c(0, max_delay_limit))+
  coord_flip() +
  theme_bw() + 
  theme(text = element_text(size = FONT_SIZE))

demo_pat_plot = ggarrange(E, F_, D, H , ncol=1, nrow=4, align = "v",
                          labels=c("A", "B", "", "C"), heights = c(1.1,0.9,0.5,0.9))

print(demo_pat_plot)
ggsave(demo_pat_plot, filename = "report/clinical/BJUI/images/demo_pat1.eps",
       device = cairo_ps, height = 6.5, width=6)
