library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

source('src/clinical_gap3/plots/dynriskPlot.R')
#do not change this order of sourcing files
load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")
source("src/clinical_gap3/fixedAndPRIASSchedule.R")
source("src/lastpaper/pers_schedule_api.R")
#this is needed to load the function which makes a plot for dynamic risk

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

FONT_SIZE = 13

schedulePlotSupp = function(object, pat_data, latest_survival_time, 
                            xbreaks, psa_breaks = NA,
                            max_follow_up, title="", M=750){
  set.seed(2019)
  
  cur_visit_time = max(pat_data$year_visit)
  
  common = ggplot() +  
    scale_x_continuous(breaks = xbreaks,
                       limits = c(-0.3, max_follow_up))+
    xlab("Follow-up time (years)")+
    geom_vline(xintercept = latest_survival_time, color=SUCCESS_COLOR) + 
    geom_vline(xintercept = cur_visit_time, linetype="dashed") +
    theme_bw() +
    theme(text = element_text(size = FONT_SIZE),
          axis.title.x = element_blank(), axis.text.x = element_blank()) 
  
  A = dynamicRiskPlot(object, pat_data, latest_survival_time, 
                      xbreaks, xbreaks, psa_breaks,
                      max_follow_up) + 
    theme(text = element_text(size = FONT_SIZE), 
          axis.title.x = element_blank(), axis.text.x = element_blank()) +
    scale_x_continuous(breaks = xbreaks,
                       limits = c(-0.3, max_follow_up)) +
    ggtitle(title)
  
  set.seed(2019)
  planned_biennial_schedule = getFixedSchedule(cur_visit_time = cur_visit_time, 
                                       latest_survival_time = latest_survival_time,
                                       min_biopsy_gap = 1, biopsy_frequency = 2, 
                                       horizon = max_follow_up)
  
  if(max(planned_biennial_schedule)<max_follow_up){
    planned_biennial_schedule = c(planned_biennial_schedule, max_follow_up)
  }
  
  planned_prias_schedule = getPRIASSchedule(object = object,
                                            patient_data = pat_data, 
                                            cur_visit_time = cur_visit_time, 
                                            latest_survival_time = latest_survival_time,
                                            M=M, horizon = max_follow_up)
  
  if(max(planned_prias_schedule)<max_follow_up){
    planned_prias_schedule = c(planned_prias_schedule, max_follow_up)
  }
  
  prias_schedule = testScheduleConsequences.mvJMbayes(object = object,
                                                      newdata = pat_data, idVar = "P_ID", 
                                                      last_test_time = latest_survival_time,
                                                      planned_prias_schedule, seed = 2019)
  biennial_schedule = testScheduleConsequences.mvJMbayes(object = object,
                                                      newdata = pat_data, idVar = "P_ID", 
                                                      last_test_time = latest_survival_time,
                                                      planned_biennial_schedule, seed = 2019)
  
  all_risk_schedules = personalizedSchedule.mvJMbayes(object = object,
                                                      newdata = pat_data, idVar = "P_ID", 
                                                      last_test_time = latest_survival_time,
                                                      fixed_grid_visits = seq(cur_visit_time, max_follow_up, 0.5),
                                                      gap = 1, seed = 2019)
  
  risk_5_perc_schedule = all_risk_schedules$all_schedules$`Risk: 0.05`
  risk_10_perc_schedule = all_risk_schedules$all_schedules$`Risk: 0.1`
  risk_15_perc_schedule = all_risk_schedules$all_schedules$`Risk: 0.15`
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
  schedule_df = rbind(schedule_df, data.frame(Schedule="Biennial", 
                                              number=1,
                                              biopsy_times=biennial_schedule$planned_test_schedule,
                                              expected_num_tests=biennial_schedule$expected_num_tests,
                                              expected_detection_delay=biennial_schedule$expected_detection_delay))
  schedule_df = rbind(schedule_df, data.frame(Schedule="PRIAS", 
                                              number=2,
                                              biopsy_times=prias_schedule$planned_test_schedule,
                                              expected_num_tests=prias_schedule$expected_num_tests,
                                              expected_detection_delay=prias_schedule$expected_detection_delay))
  schedule_df = rbind(schedule_df, data.frame(Schedule="15% Risk", 
                                              number=4,
                                              biopsy_times=risk_15_perc_schedule$planned_test_schedule,
                                              expected_num_tests=risk_15_perc_schedule$expected_num_tests,
                                              expected_detection_delay=risk_15_perc_schedule$expected_detection_delay))
  
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
  
  B =  common +
    geom_line(data = schedule_df, aes(x=biopsy_times, 
                                      y=as.numeric(schedule_df$Schedule),
                                      group=as.numeric(schedule_df$Schedule)),
              linetype='dotted')+
    geom_label(data = schedule_df, 
               aes(x=biopsy_times, y=as.numeric(schedule_df$Schedule)),
               label="B",size=POINT_SIZE + 2, fill=DANGER_COLOR, color='white') + 
    ylab("Biopsy schedule") + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size=FONT_SIZE),
          plot.margin = margin(b = 0, unit = "pt")) +
    scale_y_continuous(breaks = 1:length(levels(schedule_df$Schedule)),
                       labels = levels(schedule_df$Schedule),
                       limits = c(0.5, length(levels(schedule_df$Schedule)) + 0.5)) 
  
  C = common +
    geom_label(aes(x=c(0,latest_survival_time,cur_visit_time), y=c(0,0,0), 
                   label = c("Start AS\nGleason\ngrade 1", 
                             "Biopsy\nGleason\ngrade 1",
                             "Current\nVisit")), color='white',
               size= LABEL_SIZE, nudge_x = c(0, 0,0),
               fill=c(WARNING_COLOR, SUCCESS_COLOR, 'black')) +
    theme(text = element_text(size = FONT_SIZE),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank()) + 
    ylim(-0.25,0.25)
  
  consequences_df = schedule_df[!duplicated(schedule_df$Schedule),]
  max_delay_limit = 2
  D = ggplot() + geom_col(aes(x=rep(consequences_df$Schedule,2), 
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
  
  plot = ggarrange(A, B, C, D, ncol=1, nrow=4, align="v", labels = c("A", "B", "", "C"), 
                   common.legend = T, legend = "none", heights = c(1,1,0.35,1)) 
  
  return(plot)
}

pat1_data = prias_long_final[prias_long_final$P_ID==956 & prias_long_final$year_visit<=2,]
pat1_data$year_visit[nrow(pat1_data)] = ceiling(pat1_data$year_visit[nrow(pat1_data)])
psa_breaks = getpsaBreaks(pat1_data)

schedule_plot_pat1 = schedulePlotSupp(object = mvJoint_psa_time_scaled, 
                                      pat_data = pat1_data,
                                      latest_survival_time = 1,
                                      xbreaks = c(0, 1, 2, 3, 4, 5, 6),
                                      psa_breaks = psa_breaks,
                                      max_follow_up = 6, title = "Real Patient 1")
ggsave(schedule_plot_pat1, filename = "report/clinical/images/demo_pat1_supp.eps",
       device = cairo_ps, height = 8, width = 7)

#############
pat2_data = prias_long_final[prias_long_final$P_ID==102 & prias_long_final$year_visit<=3,]
pat2_data$year_visit[nrow(pat2_data)] = 2.5
psa_breaks = getpsaBreaks(pat2_data)

schedule_plot_pat2 = schedulePlotSupp(mvJoint_psa_time_scaled, pat2_data,
                                      latest_survival_time = 1,
                                      xbreaks = c(0, 1, 2, 3, 4,5,6),
                                      psa_breaks = psa_breaks,
                                      max_follow_up = 6, title = "Real Patient 2")

ggsave(schedule_plot_pat2, filename = "report/clinical/images/demo_pat2_supp.eps",
       device = cairo_ps, height = 8, width = 7)
##########

pat3_data = prias_long_final[prias_long_final$P_ID==101 & prias_long_final$year_visit<=4,]
pat3_data$year_visit[nrow(pat3_data)] = 4
psa_breaks = getpsaBreaks(pat3_data)

schedule_plot_pat3 = schedulePlotSupp(mvJoint_psa_time_scaled, pat3_data,
                                      latest_survival_time = 3,
                                      xbreaks = 0:6,
                                      psa_breaks = psa_breaks,
                                      max_follow_up = 6, title = "Real Patient 3")
ggsave(schedule_plot_pat3, filename = "report/clinical/images/demo_pat3_supp.eps",
       device = cairo_ps, height = 8, width = 7)

