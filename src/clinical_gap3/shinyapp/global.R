load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/gap3/PRIAS_2019/mvJoint_psa.Rdata")
print("Done")
source("prediction_only_psa.R")
source("riskSchedule.R")

#load("appdata.Rdata")
MAX_FOLLOW_UP = 10
YEAR_DIVIDER = 24 * 60 * 60 * 365
SPSS_ORIGIN_DATE = "1582-10-14"
DATE_PRINT_FORMAT = "%b %e, %Y"
DATE_PRINT_FORMAT_ABBREVIATED = "%b %y"
FONT_SIZE=15
POINT_SIZE=4
M=500

DANGER_COLOR = "#c71c22"
THEME_COLOR = "#2fa4e7"
THEME_COLOR_DARK = "#033a6f"

getHumanReadableDate = function(spss_date, abbreviated=F){
  format(as.POSIXct(spss_date, origin = SPSS_ORIGIN_DATE), 
         format = ifelse(abbreviated, DATE_PRINT_FORMAT_ABBREVIATED, DATE_PRINT_FORMAT))
}

psaObsDataGraph = function(data){
  last_visit_time = data$year_max_followup[1]
  biopsy_times = data$year_visit[!is.na(data$gleason_sum)]
  
  pred_times = seq(0, last_visit_time, length.out=50)
  
  pred_log2psaplus1 = getExpectedFutureOutcomes(mvJoint_psa, data,
                            data$latest_survival_time[1],
                            data$earliest_failure_time[1],
                            psa_predict_times = pred_times,
                            M = M)$predicted_psa
  pred_mean_psa = rowMeans(2^pred_log2psaplus1 - 1)
  
  total_x_points = 6
  x_points = seq(0, last_visit_time, length.out = total_x_points)
  x_labels = round(x_points, 1)
  x_labels[total_x_points] = paste0(x_labels[total_x_points], 
                                    "\n(Current Visit)")
  
  plot = ggplot() +
    geom_point(aes(x=data$year_visit,y=data$psa,
                   shape="Observed"), 
               size=POINT_SIZE, color=THEME_COLOR)+
    geom_line(aes(x=pred_times,y=pred_mean_psa,
              linetype="Fitted"), color=THEME_COLOR_DARK) +
    geom_segment(aes(x=biopsy_times, xend=biopsy_times,
                     y=-rep(Inf,3), yend=rep(Inf,3), linetype="Biopsy")) + 
    scale_shape_manual(name="", values = 16)+
    scale_linetype_manual(name="", values = c("solid","dashed"))+
    scale_x_continuous(breaks = x_points, labels = x_labels,
                       limits = c(0, 1.1*last_visit_time))+
    theme_bw() + 
    theme(axis.text = element_text(size = FONT_SIZE),
          axis.title = element_text(size = FONT_SIZE),
      legend.position = "bottom", legend.direction = "horizontal") + 
    xlab("Follow-up time (years)") + 
    ylab("PSA (ng/mL)")
  return(plot)
}

cumRiskGraph = function(data){
  
  last_visit_time = data$year_max_followup[1]
  if(last_visit_time >= MAX_FOLLOW_UP){
    return(NULL)
  }
  
  latest_survival_time = data$latest_survival_time[1]
  
  pred_times = seq(last_visit_time, MAX_FOLLOW_UP, length.out=50)
  
  xTicks = seq(last_visit_time, MAX_FOLLOW_UP, by=1)
  xTicks_spps_dates = data$dom_diagnosis[1] + xTicks*YEAR_DIVIDER
  xlabels = sapply(xTicks_spps_dates, getHumanReadableDate, abbreviated=T)
  xlabels[1] = paste0(getHumanReadableDate(xTicks_spps_dates[1]),"\n(Current Visit)")
  
  pred_surv_prob = getExpectedFutureOutcomes(mvJoint_psa, data,
                                             data$latest_survival_time[1],
                                             data$earliest_failure_time[1],
                                             survival_predict_times = pred_times,
                                             M = M)$predicted_surv_prob
  
  mean_risk_prob = apply(1-pred_surv_prob,1, mean)
  lower_risk_prob = apply(1-pred_surv_prob,1, quantile, probs=c(0.025))
  upper_risk_prob = apply(1-pred_surv_prob,1, quantile, probs=c(0.975))
  
  cumRiskPlot = ggplot() + 
    geom_line(aes(x=pred_times, y=mean_risk_prob), color=DANGER_COLOR) + 
    geom_ribbon(aes(x=pred_times, ymin=lower_risk_prob,
                    ymax=upper_risk_prob), fill=DANGER_COLOR, alpha=0.25) + 
    theme_bw() + 
    theme(axis.text = element_text(size = FONT_SIZE),
          axis.title = element_text(size = FONT_SIZE),
          legend.position = "bottom", legend.direction = "horizontal") + 
    scale_x_continuous(breaks = xTicks, labels=xlabels,
                       limits = c(last_visit_time, MAX_FOLLOW_UP))+
    scale_y_continuous(breaks = seq(0,1, 0.25), 
                       labels= paste0(seq(0,1, 0.25)*100, "%"),
                      limits = c(0,1)) +
    xlab("Future Dates") + 
    ylab("Cumulative Risk of Gleason Upgrade (%)")
  
  return(cumRiskPlot)
}

riskGaugeGraph = function(data, visit_time){
  
  pred_surv_prob = getExpectedFutureOutcomes(mvJoint_psa, data,
                            data$latest_survival_time[1],
                            data$earliest_failure_time[1],
                            survival_predict_times = visit_time,
                            M = M)$predicted_surv_prob
  
  mean_risk_prob = rowMeans(1-pred_surv_prob)
  risk_label = paste0(round(mean_risk_prob*100,2),"%")
  
  riskGauge = ggplot(data = NULL, 
                     aes(ymax = mean_risk_prob, ymin = 0, xmax = 2, xmin = 1, 
                     fill="Risk")) +
    geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="white", color=DANGER_COLOR) +
    geom_rect() +
    coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2.5)) + ylim(c(0,2)) +
    geom_text(aes(x = 0, y = 0, label = risk_label), 
              color=DANGER_COLOR, size=6) +
    scale_fill_manual("", values=DANGER_COLOR)+
    theme_void() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    guides(fill=FALSE) +
    guides(colour=FALSE)
  return(riskGauge)
}


