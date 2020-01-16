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
POINT_SIZE = 2
FONT_SIZE = 12
LABEL_SIZE = 2.5

set.seed(2019)
pat_data = prias_long_final[prias_long_final$P_ID==102,]
##I am perturbing the PSA of one of the patients to demo effect of rising PSA
pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] =  pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] + runif(n = 2, 0.1, 0.5)
pat_data = pat_data[pat_data$year_visit<=3.0,]

schedules = personalizedSchedule.mvJMbayes(object=mvJoint_dre_psa_2knots_quad_age,
                                           newdata = pat_data, idVar = "P_ID", last_test_time = 1.5,
                                           gap = 1, horizon = 10, seed = 2019, M = 400, cache_size = 1000)

total_schedules = length(schedules$all_schedules)
risk_thresholds = sapply(schedules$all_schedules, "[[", "threshold")
expected_delays = sapply(schedules$all_schedules, "[[", "expected_delay")
total_tests = sapply(lapply(schedules$all_schedules, "[[", "practical_test_times"), length)
euclidean_distance = sapply(schedules$all_schedules, "[[", "euclidean_distance")

min_dist_schedule_index = which.min(euclidean_distance)[1]

kappa_choice = ggplot() + 
  geom_hline(yintercept = 1.5, linetype='dashed', color=WARNING_COLOR) +
  geom_label(aes(x=6, y=1.5, label="Clinically acceptable limit for maximum time delay"), color=WARNING_COLOR, size=LABEL_SIZE) +
  geom_segment(aes(x=1,xend=total_tests[-min_dist_schedule_index], 
                   y=0,yend=expected_delays[-min_dist_schedule_index]), 
               alpha=0.125, color='gray') +
  geom_segment(aes(x=1,xend=total_tests[min_dist_schedule_index], 
                   y=0,yend=expected_delays[min_dist_schedule_index]),
               color=SUCCESS_COLOR) +
  geom_point(aes(x=total_tests[-min_dist_schedule_index], 
                 y=expected_delays[-min_dist_schedule_index]), 
             size=POINT_SIZE) +
  geom_point(aes(x=total_tests[min_dist_schedule_index], 
                 y=expected_delays[min_dist_schedule_index]), 
             size=POINT_SIZE+1, color=SUCCESS_COLOR, shape=17) +
  geom_label(aes(x=total_tests[min_dist_schedule_index], 
                 y=expected_delays[min_dist_schedule_index], 
                 label=paste0("Personalized\nSchedule\n(k = ", 
                             round(risk_thresholds[min_dist_schedule_index]*100,1), "%)")), 
             nudge_x = 1, fill=SUCCESS_COLOR, color='white', size=LABEL_SIZE)+
  geom_label(aes(x=total_tests[c(1, total_schedules)], 
                 y=expected_delays[c(1, total_schedules)], 
                 label=paste0("Personalized\nSchedule\n(k = ", 
                             round(risk_thresholds[c(1, total_schedules)]*100,1), "%)")), 
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