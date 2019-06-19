load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
source("prediction_only_psa.R")
source("riskSchedule.R")
source("schedules.R")

#load("appdata.Rdata")
MAX_FOLLOW_UP = 10
YEAR_DIVIDER = 24 * 60 * 60 * 365
SPSS_ORIGIN_DATE = "1582-10-14"
DATE_PRINT_FORMAT = "%b %e, %Y"
DATE_PRINT_FORMAT_ABBREVIATED = "%b %Y"
FONT_SIZE=15
POINT_SIZE=4
M=500
CACHE_SIZE=200

DANGER_COLOR = "#c71c22"
THEME_COLOR = "#2fa4e7"
THEME_COLOR_DARK = "#033a6f"
SUCCESS_COLOR = "#73A839"
WARNING_COLOR="orange"

getHumanReadableDate = function(spss_date, abbreviated=F){
  format(as.POSIXct(spss_date, origin = SPSS_ORIGIN_DATE), 
         format = ifelse(abbreviated, DATE_PRINT_FORMAT_ABBREVIATED, DATE_PRINT_FORMAT))
}

psaObsDataGraph = function(data, pred_times, pred_log2psaplus1){
  last_visit_time = data$year_max_followup[1]
  biopsy_times = data$year_visit[!is.na(data$gleason_sum)]
  
  pred_mean_psa = rowMeans(2^pred_log2psaplus1 - 1)
  
  total_x_points = 6
  x_points = seq(0, last_visit_time, length.out = total_x_points)
  x_labels = round(x_points, 1)
  x_labels[total_x_points] = paste0(x_labels[total_x_points], 
                                    "\n(Latest Visit)")
  
  plot = ggplot() +
    geom_ribbon(aes(x=c(0, data$latest_survival_time[1]),
                    ymin=-c(Inf,Inf), ymax=c(Inf,Inf), fill="Region with\nGleason ≤ 6"),
                alpha=0.25) +
    geom_point(aes(x=data$year_visit,y=data$psa,
                   shape="Observed PSA"), 
               size=POINT_SIZE, color=THEME_COLOR)+
    geom_line(aes(x=pred_times,y=pred_mean_psa,
                  linetype="Fitted PSA"), color=THEME_COLOR_DARK) +
    geom_segment(aes(x=biopsy_times, xend=biopsy_times,
                     y=-rep(Inf,length(biopsy_times)), 
                     yend=rep(Inf,length(biopsy_times)), linetype="Biopsy")) + 
    scale_fill_manual(name="", values = SUCCESS_COLOR) + 
    scale_shape_manual(name="", values = 16)+
    scale_linetype_manual(name="", values = c("solid","dashed"))+
    scale_x_continuous(breaks = x_points, labels = x_labels,
                       limits = c(0, 1.1*last_visit_time))+
    theme_bw() + 
    theme(axis.text = element_text(size = FONT_SIZE),
          axis.title = element_text(size = FONT_SIZE),
          legend.text = element_text(size = FONT_SIZE),
          legend.position = "bottom", legend.direction = "horizontal") + 
    xlab("Follow-up time (years)") + 
    ylab("PSA (ng/mL)")
  return(plot)
}

biopsyScheduleGraph = function(schedule_df, 
                               selected_schedule_ids,
                               start_dom,
                               current_visit_time,
                               latest_survival_time){
  
  schedule_df = schedule_df[schedule_df$schedule_id %in% selected_schedule_ids,]
  
  xTicks = seq(0, MAX_FOLLOW_UP, length.out = 6)
  xTicks_spps_dates = start_dom + xTicks*YEAR_DIVIDER
  xlabels = sapply(xTicks_spps_dates, getHumanReadableDate, abbreviated=T)
  xlabels[1] = paste0(getHumanReadableDate(xTicks_spps_dates[1]),"\n(Start AS)")
  
  ybreaks = unique(schedule_df$schedule_id)
  ylabels = unique(schedule_df$Schedule)
  
  pp = ggplot() +
    geom_ribbon(aes(x=c(0, latest_survival_time),
                    ymin=-c(Inf,Inf), ymax=c(Inf,Inf), fill="Region with Gleason ≤ 6"),
                alpha=0.25) +
    geom_vline(xintercept = latest_survival_time) +
    geom_vline(xintercept = current_visit_time, linetype='dashed') +
    geom_label(data=schedule_df,
               aes(x=biopsy_times, y=schedule_id, label='B', group=schedule_id), color=THEME_COLOR,
               size=7) +
    geom_line(data=schedule_df,
              aes(x=biopsy_times, y=schedule_id, group=schedule_id), color=THEME_COLOR,
              linetype='dotted') +
    scale_fill_manual(name="", values = SUCCESS_COLOR) + 
    theme_bw() + xlab("Proposed Dates of Biopsy") +
    scale_x_continuous(breaks = xTicks, labels=xlabels,
                       limits = c(0, MAX_FOLLOW_UP))+
    scale_y_continuous(breaks = ybreaks, labels = ylabels, 
                       limits = c(ybreaks[1]-0.5, tail(ybreaks,1) + 0.5)) +
    ylab("Biopsy Schedule") + 
    theme(text = element_text(size = FONT_SIZE),
          legend.position = "bottom", legend.direction = "horizontal")
  
  return(pp)
}

biopsyDelayGraph = function(schedule_df, 
                            selected_schedule_ids){
  schedule_df = schedule_df[schedule_df$schedule_id %in% selected_schedule_ids,]
  ybreaks = seq(0,ceiling(max(schedule_df$expected_delay)), length.out = 4)
  ylabels = round(ybreaks,1)
  
  pp = ggplot(data=schedule_df) + geom_bar(aes(x=schedule, y=expected_delay),
                          stat='identity', width=0.5) +  
    ylab("Expected delay (months)\n in detecting Gleason \u2265 7") + 
    xlab("Biopsy schedule") +
    scale_y_continuous(breaks = ybreaks,
                       labels = ylabels,
                       limits = c(ybreaks[1], tail(ybreaks,1)))+
    coord_flip() +
    theme_bw() + 
    theme(text = element_text(size = FONT_SIZE),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  return(pp)
}

biopsyTotalGraph = function(schedule_df, 
                            selected_schedule_ids){  
  schedule_df = schedule_df[schedule_df$schedule_id %in% selected_schedule_ids,]
  
  pp = ggplot(data=schedule_df) + geom_bar(aes(x=schedule, y=total_biopsies),
                          stat='identity', width=0.5) +  
    ylab("Total biopsies scheduled") + 
    xlab("Biopsy schedule") + 
    theme_bw() + 
    theme(text = element_text(size = FONT_SIZE))+
    coord_flip()
  return(pp)
}

riskGaugeGraph = function(mean_risk_prob, date=""){
  
  gauge_color = colorRamp(c(SUCCESS_COLOR, WARNING_COLOR, DANGER_COLOR))(mean_risk_prob)
  gauge_color = rgb(gauge_color[1], gauge_color[2], gauge_color[3], maxColorValue = 255)
  
  risk_label = paste0("\n\n\nCumulative-risk: ", round(mean_risk_prob*100,2),"%\n on ",date)
  
  riskGauge = ggplot(data = NULL, 
                     aes(ymax = mean_risk_prob, ymin = 0, xmax = 2, xmin = 1, 
                         fill="Risk")) +
    geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="white", color=gauge_color) +
    geom_rect() +
    coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2.5)) + ylim(c(0,2)) +
    geom_text(aes(x = 0, y = 0, label = risk_label), 
              color=gauge_color, size=6) +
    geom_segment(aes(x=0, xend=1, y=mean_risk_prob, yend=mean_risk_prob), 
                 color=gauge_color,
                 arrow = arrow(length = unit(0.4,"cm")))+
    geom_point(aes(x=0, y=0), color=gauge_color, size=3)+
    scale_fill_manual("", values=gauge_color)+
    theme_void() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size=15, color=gauge_color)) +
    guides(fill=FALSE) +
    guides(colour=FALSE)
  return(riskGauge)
}


