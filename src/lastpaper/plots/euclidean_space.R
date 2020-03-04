library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_quad_age.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/lastpaper/pers_schedule_api.R")

SUCCESS_COLOR = 'forestgreen'
DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10
POINT_SIZE = 4
FONT_SIZE = 12
LABEL_SIZE = 2.5

set.seed(2019)
pat_data = prias_long_final[prias_long_final$P_ID==102,]
##I am perturbing the PSA of one of the patients to demo effect of rising PSA
pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] =  pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] + runif(n = 2, 0.1, 0.5)
pat_data = pat_data[pat_data$year_visit<=3.0,]
pat_data$year_visit[nrow(pat_data)] = 2.5
cur_visit_time = pat_data$year_visit[nrow(pat_data)]

schedules = personalizedSchedule.mvJMbayes(object=mvJoint_dre_psa_2knots_quad_age,
                                           newdata = pat_data, idVar = "P_ID", last_test_time = 1.5,
                                           gap = 1, fixed_grid_visits = seq(cur_visit_time, to = MAX_FOLLOW_UP, by = 0.5), 
                                           seed = 2019)

total_schedules = length(schedules$all_schedules)
risk_thresholds = sapply(schedules$all_schedules, "[[", "cumulative_risk_threshold")
expected_delays = sapply(schedules$all_schedules, "[[", "expected_detection_delay")
expected_total_tests = sapply(schedules$all_schedules, "[[", "expected_num_tests")
euclidean_distance = sapply(schedules$all_schedules, "[[", "euclidean_distance")

ceiling_expected_delay = ceiling(max(expected_delays))
ceiling_expected_total_tests = ceiling(max(expected_total_tests))

min_dist_schedule_index = which.min(euclidean_distance)[1]

delay_threshold = 1.5

kappa_choice = ggplot() + 
  geom_hline(yintercept = delay_threshold, linetype='dashed', color=WARNING_COLOR) +
  geom_label(aes(x=3.5, y=delay_threshold, label="Clinically acceptable limit for maximum time delay (example)"), color=WARNING_COLOR, size=LABEL_SIZE) +
  geom_segment(aes(x=1,xend=expected_total_tests[-min_dist_schedule_index], 
                   y=0,yend=expected_delays[-min_dist_schedule_index]), 
               alpha=0.175, color='gray') +
  geom_segment(aes(x=1,xend=expected_total_tests[min_dist_schedule_index], 
                   y=0,yend=expected_delays[min_dist_schedule_index]),
               color=SUCCESS_COLOR) +
  geom_point(aes(x=expected_total_tests[-min_dist_schedule_index], 
                 y=expected_delays[-min_dist_schedule_index]), 
             size=POINT_SIZE) +
  geom_point(aes(x=expected_total_tests[min_dist_schedule_index], 
                 y=expected_delays[min_dist_schedule_index]), 
             size=POINT_SIZE+2, color=SUCCESS_COLOR, shape=17) +
  geom_label(aes(x=expected_total_tests[min_dist_schedule_index], 
                 y=expected_delays[min_dist_schedule_index], 
                 label=paste0("Personalized\nSchedule\n\u03BA*(v) = ", 
                              round(risk_thresholds[min_dist_schedule_index]*100,1), "%")), 
             nudge_x = 0, nudge_y = -0.5, fill=SUCCESS_COLOR, color='white', size=LABEL_SIZE)+
  geom_label(aes(x=expected_total_tests[c(1, total_schedules)], 
                 y=expected_delays[c(1, total_schedules)], 
                 label=paste0("Personalized\nSchedule\n\u03BA = ", 
                              round(risk_thresholds[c(1, total_schedules)]*100,1), "%")), 
             nudge_x = c(0.4, 0.4), nudge_y=c(0.2, 0), fill='black', color='white', size=LABEL_SIZE)+
  geom_point(aes(x=1, y=0), shape=15, size=POINT_SIZE + 1, 
             color=THEME_COLOR) +
  geom_label(aes(x=1,y=0, label="Ideal Schedule"), 
             nudge_y = -0.25,
             fill=THEME_COLOR, color='white',
             size=LABEL_SIZE)+
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE)) +
  scale_x_continuous(breaks=1:ceiling_expected_total_tests, 
                     limits = c(0.5,ceiling_expected_total_tests)) +
  scale_y_continuous(breaks=seq(0, ceiling_expected_delay, 1), 
                     limits = c(-0.5,ceiling_expected_delay)) +
  xlab("Expected number of tests") +
  ylab("Expected time delay in detecting progression (e.g., months,years)")

print(kappa_choice)
ggsave(kappa_choice, filename = "report/lastpaper/images/kappa_choice_102.eps",
       device = cairo_ps,  height = 6, width=6)
