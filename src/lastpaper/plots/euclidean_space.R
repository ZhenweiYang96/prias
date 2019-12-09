library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")
source("src/lastpaper/minDistThreshold.R")
source("src/lastpaper/scheduleCreatorCacheBased.R")

SUCCESS_COLOR = 'forestgreen'
DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10
POINT_SIZE = 2
FONT_SIZE = 12
LABEL_SIZE = 2.5

set.seed(2019)
pat_data = prias_long_final[prias_long_final$P_ID==102,]
##I am perturbing the PSA of one of the patients to demo effect of rising PSA
pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] =  pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] + runif(n = 2, 0.1, 0.5)

schedules = minDistScheduleDecision(object = mvJoint_psa_time_scaled, 
                             patient_data = pat_data[pat_data$year_visit<=3.0,],
                             cur_visit_time = 2.5890411,
                             latest_survival_time = 1.5)
min_dist_schedule_index = which.min(schedules$dist)
total_schedules = length(schedules$dist)

kappa_choice = ggplot() + 
  geom_hline(yintercept = 1.5, linetype='dashed', color=WARNING_COLOR) +
  geom_label(aes(x=6, y=1.5, label="Clinically acceptable threshold for maximum time delay"), color=WARNING_COLOR, size=LABEL_SIZE) +
  geom_segment(aes(x=1,xend=schedules$total_biopsies[-min_dist_schedule_index], 
                   y=0,yend=schedules$expected_delays[-min_dist_schedule_index]), 
               alpha=0.125, color='gray') +
  geom_segment(aes(x=1,xend=schedules$total_biopsies[min_dist_schedule_index], 
                   y=0,yend=schedules$expected_delays[min_dist_schedule_index]),
               color=SUCCESS_COLOR) +
  geom_point(aes(x=schedules$total_biopsies[-min_dist_schedule_index], 
                 y=schedules$expected_delays[-min_dist_schedule_index]), 
             size=POINT_SIZE) +
  geom_point(aes(x=schedules$total_biopsies[min_dist_schedule_index], 
                 y=schedules$expected_delays[min_dist_schedule_index]), 
             size=POINT_SIZE+1, color=SUCCESS_COLOR, shape=17) +
  geom_label(aes(x=schedules$total_biopsies[min_dist_schedule_index], 
                 y=schedules$expected_delays[min_dist_schedule_index], 
                 label=paste0("Personalized\nSchedule\n(k = ", 
                             round(schedules$risk_thresholds[min_dist_schedule_index]*100,1), "%)")), 
             nudge_x = 1, fill=SUCCESS_COLOR, color='white', size=LABEL_SIZE)+
  geom_label(aes(x=schedules$total_biopsies[c(1, total_schedules)], 
                 y=schedules$expected_delays[c(1, total_schedules)], 
                 label=paste0("Personalized\nSchedule\n(k = ", 
                             round(schedules$risk_thresholds[c(1, total_schedules)]*100,1), "%)")), 
             nudge_x = c(0, 1), nudge_y=c(0.2, 0), fill='black', color='white', size=LABEL_SIZE)+
  geom_point(aes(x=1, y=0), shape=15, size=POINT_SIZE + 1, 
             color=THEME_COLOR) +
  geom_label(aes(x=1,y=0, label="Ideal Schedule"), 
             nudge_y = -0.1,
             fill=THEME_COLOR, color='white',
             size=LABEL_SIZE)+
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE)) +
  scale_x_continuous(breaks=1:10, limits = c(0.5,9)) +
  scale_y_continuous(breaks=seq(0, 2.5, 0.5), limits = c(-0.2,2.25)) +
  xlab("Number of biopsies") +
  ylab("Expected time delay in detecting progression")

print(kappa_choice)
ggsave(kappa_choice, filename = "report/lastpaper/images/kappa_choice_102.eps",
       device = cairo_ps,  height = 5, width=5)
