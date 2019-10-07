library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")
#this is needed to load the function which makes a plot for dynamic risk
source('src/clinical_gap3/plots/dynriskPlot.R')
source('src/clinical_gap3/scheduleCreator.R')

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
  M=750
  schedule_5perc = getRiskBasedSchedule(object = mvJoint_psa_time_scaled, patient_data = pat_data,
                                        cur_visit_time = cur_visit_time, 
                                        latest_survival_time = latest_survival_time,
                                        risk_threshold = 0.05, min_biopsy_gap = 1, M = M,
                                        horizon = MAX_FOLLOW_UP)
  schedule_10perc = getRiskBasedSchedule(object = mvJoint_psa_time_scaled, patient_data = pat_data,
                                         cur_visit_time = cur_visit_time, 
                                         latest_survival_time = latest_survival_time,
                                         risk_threshold = 0.1, min_biopsy_gap = 1, M = M,
                                         horizon = MAX_FOLLOW_UP)
  schedule_15perc = getRiskBasedSchedule(object = mvJoint_psa_time_scaled, patient_data = pat_data,
                                         cur_visit_time = cur_visit_time, 
                                         latest_survival_time = latest_survival_time,
                                         risk_threshold = 0.15, min_biopsy_gap = 1, M = M,
                                         horizon = MAX_FOLLOW_UP)
  schedule_prias = getPRIASSchedule(object = mvJoint_psa_time_scaled, patient_data = pat_data,
                                    cur_visit_time = cur_visit_time, latest_survival_time = latest_survival_time,
                                    min_biopsy_gap = 1, M = M,horizon = MAX_FOLLOW_UP)
  schedule_annual = getFixedSchedule(cur_visit_time = cur_visit_time, 
                                     latest_survival_time = latest_survival_time,
                                     min_biopsy_gap = 1, biopsy_frequency = 1, 
                                     horizon = MAX_FOLLOW_UP)
  schedule_biennial = getFixedSchedule(cur_visit_time = cur_visit_time, 
                                       latest_survival_time = latest_survival_time,
                                       min_biopsy_gap = 1, biopsy_frequency = 2, 
                                       horizon = MAX_FOLLOW_UP)
  
  schedules = c("5% Risk", "10% Risk", "15% Risk","PRIAS", "Yearly", "Biennial")
  
  consequences_list = lapply(list(schedule_5perc, schedule_10perc,schedule_15perc,
                                  schedule_prias, schedule_annual,schedule_biennial),
                             FUN = function(x){
                               set.seed(2019); 
                               getConsequences(object=mvJoint_psa_time_scaled, 
                                               patient_data=pat_data,
                                               cur_visit_time = cur_visit_time, 
                                               proposed_biopsy_times = x,
                                               latest_survival_time = latest_survival_time, M=M,
                                               horizon=MAX_FOLLOW_UP)})
  
  expected_delays = sapply(consequences_list, "[[", "expected_delay")
  practical_biopsy_times = lapply(consequences_list, "[[", "practical_biopsy_times")
  total_biopsies = sapply(practical_biopsy_times, length)
  
  plotDf = data.frame(Schedule=rep(schedules, total_biopsies), 
                      biopsy_times=unlist(practical_biopsy_times))  
  B = common +
    geom_line(data = plotDf, aes(x=biopsy_times, 
                                 y=as.numeric(plotDf$Schedule),
                                 group=as.numeric(plotDf$Schedule)),
              linetype='dotted')+
    geom_label(data = plotDf, 
               aes(x=biopsy_times, y=as.numeric(plotDf$Schedule)),
               label="B",size=POINT_SIZE + 2, fill=DANGER_COLOR, color='white') + 
    ylab("Biopsy schedule") + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size=FONT_SIZE),
          plot.margin = margin(b = 0, unit = "pt")) +
    scale_y_continuous(breaks = 1:length(levels(plotDf$Schedule)),
                       labels = levels(plotDf$Schedule),
                       limits = c(0.5, length(levels(plotDf$Schedule)) + 0.5)) 
  
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
  
  D = ggplot() + geom_col(aes(x=rep(schedules,2), 
                              y=c(12*expected_delays, 12 - 12*expected_delays)),
                          color='black', fill=c(rep(c('darkgrey','white'), length(schedules))),
                          width=0.5)+
    ylab("Expected delay (months) in detecting reclassification") + 
    xlab("Biopsy schedule") +
    scale_y_continuous(breaks = seq(0,12, length.out = 7),
                       labels = round(seq(0,12, length.out = 7),1),
                       limits=c(0,12))+
    coord_flip() +
    theme_bw() + 
    theme(text = element_text(size = FONT_SIZE))
  
  plot = ggarrange(A, B, C, D, ncol=1, nrow=4, align="v", labels = c("A", "B", "", "C"), 
                   common.legend = T, legend = "none", heights = c(1,1,0.35,1)) 
  
  return(plot)
}

pat1_data = prias_long_final[prias_long_final$P_ID==956 & prias_long_final$year_visit<=2,]
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
psa_breaks = getpsaBreaks(pat3_data)

schedule_plot_pat3 = schedulePlotSupp(mvJoint_psa_time_scaled, pat3_data,
                                      latest_survival_time = 3,
                                      xbreaks = 0:6,
                                      psa_breaks = psa_breaks,
                                      max_follow_up = 6, title = "Real Patient 3")
ggsave(schedule_plot_pat3, filename = "report/clinical/images/demo_pat3_supp.eps",
       device = cairo_ps, height = 8, width = 7)

