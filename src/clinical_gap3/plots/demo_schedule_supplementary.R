source('src/clinical_gap3/plots/dynriskPlot.R')
source('src/clinical_gap3/compareSchedules.R')

FONT_SIZE = 14

schedulePlotSupp = function(object, pat_df, latest_survival_time=NA, 
                            xbreaks, xlabs, psa_breaks = NA,
                            max_follow_up, title=""){
  set.seed(2019)
  
  if(is.na(latest_survival_time)){
    latest_survival_time = max(pat_df$year_visit[!is.na(pat_df$gleason_sum)])
  }
  
  max_psa_time = max(pat_df$year_visit)   
  
  A = dynamicRiskPlot(mvJoint_psa_time_scaled, 
                      pat_df, latest_survival_time, 
                      xbreaks, xlabs, psa_breaks,
                      max_follow_up) + theme(axis.title.x = element_blank())
  A = A + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank()) +
    ggtitle(title)
  
  common = ggplot() +  
    scale_x_continuous(breaks = xbreaks, labels=xlabs,
                       limits = c(0, max_follow_up), 
                       minor_breaks = seq(0, max_follow_up, 1)) +
    xlab("Follow-up time (years)")+
    geom_ribbon(aes(x=c(0,latest_survival_time), ymin=-c(Inf,Inf),
                    ymax=c(Inf, Inf)), alpha=0.15, fill=SUCCESS_COLOR) + 
    geom_vline(xintercept = latest_survival_time) + 
    geom_vline(xintercept = max_psa_time, linetype="dashed") +
    theme_bw() +
    theme(text = element_text(size = FONT_SIZE)) 
  
  biopsy_schedules = compareSchedules(pat_df,cur_visit_time = max_psa_time, 
                                      latest_survival_time = latest_survival_time)
  
  schedules = c("5% Risk", "10% Risk", "15% Risk", 
                "PRIAS", "Yearly", "Every 2 Years")
  expected_delays = sapply(biopsy_schedules$schedules, "[[", "expected_delay")
  total_biopsies = sapply(biopsy_schedules$schedules, "[[", "total_biopsies")
  
  biopsy_times = do.call('c', lapply(biopsy_schedules$schedules, "[[", "biopsy_times"))
  plotDf = data.frame(Schedule=rep(schedules, total_biopsies), 
                      biopsy_times)
  
  B = common +
    geom_line(data = plotDf, aes(x=biopsy_times, 
                                 y=as.numeric(plotDf$Schedule),
                                 group=as.numeric(plotDf$Schedule)),
              linetype='dotted')+
    geom_label(data = plotDf, 
               aes(x=biopsy_times, y=as.numeric(plotDf$Schedule)),
               label="B",size=POINT_SIZE + 2) + 
    ylab("Biopsy schedule") + 
    scale_y_continuous(breaks = 1:length(levels(plotDf$Schedule)),
                       labels = levels(plotDf$Schedule),
                       limits = c(0.5, length(levels(plotDf$Schedule)) + 0.5)) 
  
  
  D = ggplot() + geom_bar(aes(x=schedules, y=12*expected_delays),
                          stat='identity', width=0.5) +  
    ylab("Expected delay (months)\n in detecting Gleason \u2265 7") + 
    xlab("Biopsy schedule") +
    scale_y_continuous(breaks = seq(0,12*ceiling(max(expected_delays)), length.out = 4),
                       labels = round(seq(0,12*ceiling(max(expected_delays)), length.out = 4),1),
                       limits=c(0,12*ceiling(max(expected_delays))))+
    coord_flip() +
    theme_bw() + 
    theme(text = element_text(size = FONT_SIZE))
  
  E = ggplot() + geom_bar(aes(x=schedules, y=total_biopsies),
                          stat='identity', width=0.5) +  
    ylab("Total biopsies scheduled") + 
    xlab("Biopsy schedule") + 
    theme_bw() + 
    theme(text = element_text(size = FONT_SIZE),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    coord_flip()
  
  upper_plot = ggarrange(A, B, ncol=1, nrow=2, align = "v", labels = "AUTO",
                         heights = c(1,1.35), legend = "bottom", common.legend = T)
  lower_plot = ggarrange(D, E, ncol=2, nrow=1, 
                         align="h", widths = c(1.25,1))
  
  return(ggarrange(upper_plot, lower_plot, ncol=1, nrow=2, heights = c(2.5,1)))
}

pat1_data = prias_long_final[prias_long_final$P_ID==956 & prias_long_final$year_visit<=2,]
psa_breaks = getpsaBreaks(pat1_data)

schedule_plot_pat1 = schedulePlotSupp(mvJoint_psa_time_scaled, pat1_data,
                                      latest_survival_time = 0,
                                      xbreaks = c(0, 1.9780822, 5, 7.5, 10),
                                      xlabs = c("0\n(Latest\nbiopsy)", "2\n (Current\nvisit)",
                                                "5","7.5", "10"),
                                      psa_breaks = psa_breaks,
                                      max_follow_up = 10, title = "Real Patient 1") + theme(axis.title.x = element_blank())
ggsave(schedule_plot_pat1, filename = "report/clinical/images/demo_pat1_supp.eps",
       device = cairo_ps, height = 8)

#############
pat2_data = prias_long_final[prias_long_final$P_ID==102 & prias_long_final$year_visit<=3,]
psa_breaks = getpsaBreaks(pat2_data)

schedule_plot_pat2 = schedulePlotSupp(mvJoint_psa_time_scaled, pat2_data,
                                      latest_survival_time = 1,
                                      xbreaks = c(0, 1, 2.5, 7.5, 10),
                                      xlabs = c("0", "1\n(Latest\nbiopsy)",
                                                "2.5\n (Current\nvisit)", "7.5", "10"),
                                      psa_breaks = psa_breaks,
                                      max_follow_up = 10, title = "Real Patient 2") + theme(axis.title.x = element_blank())

ggsave(schedule_plot_pat2, filename = "report/clinical/images/demo_pat2_supp.eps",
       device = cairo_ps, height = 8)
##########

pat3_data = prias_long_final[prias_long_final$P_ID==101 & prias_long_final$year_visit<=5,]
psa_breaks = getpsaBreaks(pat3_data)

schedule_plot_pat3 = schedulePlotSupp(mvJoint_psa_time_scaled, pat3_data,
                                      latest_survival_time = 3,
                                      xbreaks = c(0, 3, 4.9178082, 7.5, 10),
                                      xlabs = c("0", "3\n(Latest\nbiopsy)",
                                                "5\n (Current\nvisit)", "7.5", "10"),
                                      psa_breaks = psa_breaks,
                                      max_follow_up = 10, title = "Real Patient 3") + theme(axis.title.x = element_blank())
ggsave(schedule_plot_pat3, filename = "report/clinical/images/demo_pat3_supp.eps",
       device = cairo_ps, height = 8)

##########

pat4_data = prias_long_final[prias_long_final$P_ID==156 & prias_long_final$year_visit<=7,]
psa_breaks = getpsaBreaks(pat4_data)

schedule_plot_pat4 = schedulePlotSupp(mvJoint_psa_time_scaled, pat4_data,
                                      latest_survival_time = 5.5,
                                      xbreaks = c(0, 3, 5.5, 6.5315068, 10),
                                      xlabs = c("0", "3","5.5\n(Latest\nbiopsy)",
                                                "6.5\n    (Current\n   visit)", "10"),
                                      psa_breaks = psa_breaks,
                                      max_follow_up = 10, title = "Real Patient 4") + theme(axis.title.x = element_blank())
ggsave(schedule_plot_pat4, filename = "report/clinical/images/demo_pat4_supp.eps",
       device = cairo_ps, height = 8)

##########

