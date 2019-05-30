load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/gap3/PRIAS_2019/mvJoint_psa.Rdata")
print("Done")
source("prediction_only_psa.R")

#load("appdata.Rdata")
MAX_FOLLOW_UP = 13.5
SPSS_ORIGIN_DATE = "1582-10-14"
DATE_PRINT_FORMAT = "%b %e, %Y"
FONT_SIZE=15
POINT_SIZE=4
M=500

DANGER_COLOR = "#c71c22"
THEME_COLOR = "#2fa4e7"
THEME_COLOR_DARK = "#033a6f"

getHumanReadableDate = function(spss_date){
  format(as.POSIXct(spss_date, origin = SPSS_ORIGIN_DATE), format = DATE_PRINT_FORMAT)
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


