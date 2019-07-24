load("mvJoint_psa_time_scaled_light.Rdata")
load("demo_pat_list.Rdata")
source("prediction_only_psa.R")
source("scheduleCreator.R")

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
SCHEDULES = c("5% Risk", "10% Risk", "15% Risk",
             "Yearly", "Every 2 Years", "PRIAS")
DELAY_GAUGE_MAX = 24

DANGER_COLOR = "red"
THEME_COLOR = "#2fa4e7"
THEME_COLOR_DARK = "#033a6f"
SUCCESS_COLOR = "forestgreen"
WARNING_COLOR="orange"

EXAMPLE_DF = data.frame(age=62.3,
                       start_date = "21-02-2016",
                       visit_date = c("21-02-2016", "20-08-2016", "15-02-2017", "19-08-2017", 
                                      "21-02-2018", "13-08-2018", "26-02-2019", "23-08-2019"),
                       psa = c(5.7, NA, 12, 8.5, 15, NA, 25, 20.3),
                       gleason_sum = c(6,NA, 6, NA, NA, 6, NA, NA))

getHumanReadableDate = function(spss_date, abbreviated=F){
  format(as.POSIXct(spss_date, origin = SPSS_ORIGIN_DATE), 
         format = ifelse(abbreviated, DATE_PRINT_FORMAT_ABBREVIATED, DATE_PRINT_FORMAT))
}

psaObsDataGraph = function(data, dom_diagnosis, current_visit_time, 
                           latest_survival_time,
                           pred_times, pred_log2psaplus1){
  biopsy_times = data$year_visit[!is.na(data$gleason_sum)]
  
  pred_mean_psa = rowMeans(2^pred_log2psaplus1 - 1)
  
  xTicks = seq(0, current_visit_time, length.out = 5)
  if(current_visit_time != latest_survival_time & latest_survival_time!=0){
    nearest_tick_index = c(2:4)[which.min(abs(xTicks[2:4] - latest_survival_time))]
    xTicks[nearest_tick_index] = latest_survival_time
    if(nearest_tick_index==2){
      xTicks[3:4] = seq(latest_survival_time, current_visit_time, length.out = 4)[2:3]
    }else if(nearest_tick_index==4){
      xTicks[2:3] = seq(0, latest_survival_time, length.out = 4)[2:3]
    }else{
      xTicks[2] = latest_survival_time/2
      xTicks[4] = (latest_survival_time + current_visit_time)/2
    }
  }
  
  latest_biopsy_label = NULL
  if(latest_survival_time==0){
    startAS_label = "Start AS\n& Latest\nbiopsy"
    current_visit_label = "Current\nvisit"
  }else{
    startAS_label = "Start AS"
    
    if(current_visit_time == latest_survival_time){
      current_visit_label = "Current visit\n& Latest biopsy"
    }else{
      current_visit_label = "Current\nvisit"
      latest_biopsy_label = "Latest\nbiopsy"
    }
  }
  
  xlabels = as.character(round(xTicks,1))
  xTicks_spps_dates = sapply(dom_diagnosis + xTicks*YEAR_DIVIDER, getHumanReadableDate, abbreviated=T)
  xlabels = sapply(1:length(xlabels), function(i){paste0(xlabels[i],"\n(",xTicks_spps_dates[i], ")")})
  
  psa_min = min(data$psa, na.rm = T)
  psa_max = max(data$psa, na.rm = T)

  psa_breaks = seq(psa_min, psa_max, length.out = 4)
  
  plot = ggplot() +
    geom_ribbon(aes(x=c(0, latest_survival_time),
                    ymin=-c(Inf,Inf), ymax=c(Inf,Inf), fill="Gleason \u2264 6"),
                alpha=0.25) +
    geom_point(aes(x=data$year_visit,y=data$psa,
                   shape="Observed PSA"), 
               size=POINT_SIZE, color=THEME_COLOR)+
    geom_line(aes(x=pred_times,y=pred_mean_psa,
                  linetype="Fitted PSA"), color=THEME_COLOR_DARK, size=0.75) +
    geom_segment(aes(x=biopsy_times, xend=biopsy_times,
                     y=-rep(Inf,length(biopsy_times)), 
                     yend=rep(Inf,length(biopsy_times)), linetype="Biopsy")) + 
    geom_vline(xintercept = current_visit_time, linetype='dashed') +
    geom_label(aes(x=c(0, current_visit_time), 
                   y=rep(psa_max * 1.1,2),
                   label=c(startAS_label, current_visit_label)), 
               fill=WARNING_COLOR, color='white', size=6, label.r = unit(0.25,units = "lines")) +
    scale_fill_manual(name="", values = SUCCESS_COLOR) + 
    scale_shape_manual(name="", values = 16)+
    scale_linetype_manual(name="", values = c("solid","dotted"))+
    scale_x_continuous(breaks = xTicks, labels = xlabels,
                       limits = c(-0.2, 1.2*current_visit_time))+
    scale_y_continuous(breaks = psa_breaks, labels=round(psa_breaks,1),
                       limits = c(psa_min, psa_max *1.4)) +
    theme_bw() + 
    theme(text = element_text(size = FONT_SIZE),
          axis.text.y = element_text(size = FONT_SIZE),
          axis.text.x = element_text(size = FONT_SIZE,
                                     angle = 30, hjust = 1),
          legend.position = "bottom", legend.direction = "horizontal")+
    xlab("Follow-up time (years)") + 
    ylab("PSA (ng/mL)")
  
  if(!is.null(latest_biopsy_label)){
    plot = plot + geom_label(aes(x=latest_survival_time, y=psa_max*1.3, 
                             label=latest_biopsy_label), 
                         fill=WARNING_COLOR, color='white', size=6, label.r = unit(0.25,units = "lines"))
  }
  
  return(plot)
}

biopsyScheduleGraph = function(schedule_df, 
                               selected_schedule_ids,
                               start_dom,
                               current_visit_time,
                               latest_survival_time,
                               prev_biopsies){
  
  schedule_df = schedule_df[schedule_df$schedule_id %in% selected_schedule_ids,]
  
  if(current_visit_time == latest_survival_time){
    if(MAX_FOLLOW_UP-current_visit_time >= 1){
      xTicks = c(0, seq(current_visit_time, MAX_FOLLOW_UP, length.out = 5))
    }else{
      xTicks = c(seq(0, current_visit_time, length.out = 5), MAX_FOLLOW_UP)
    }
  }else{
    if(MAX_FOLLOW_UP-current_visit_time >= 1){
      xTicks = c(0, latest_survival_time, seq(current_visit_time, MAX_FOLLOW_UP, length.out = 4))
    }else{
      xTicks = c(0,latest_survival_time/2, latest_survival_time, (latest_survival_time+current_visit_time)/2, current_visit_time, MAX_FOLLOW_UP)
    }
  }
  
  xlabels = as.character(round(xTicks,1))
  xTicks_spps_dates = sapply(start_dom + xTicks*YEAR_DIVIDER, getHumanReadableDate, abbreviated=T)
  xlabels = sapply(1:length(xlabels), function(i){paste0(xlabels[i],"\n(",xTicks_spps_dates[i], ")")})
  
  latest_biopsy_label = NULL
  if(latest_survival_time==0){
    startAS_label = "Start AS\n& Latest\nbiopsy"
    current_visit_label = "Current\nvisit"
  }else{
    startAS_label = "Start AS"
    
    if(current_visit_time == latest_survival_time){
      current_visit_label = "Current visit\n& Latest biopsy"
    }else{
      current_visit_label = "Current\nvisit"
      latest_biopsy_label = "Latest\nbiopsy"
    }
  }
  
  ybreaks = unique(schedule_df$schedule_id)
  ylabels = unique(schedule_df$Schedule)
  
  pp = ggplot() +
    geom_ribbon(aes(x=c(0, latest_survival_time),
                    ymin=-c(Inf,Inf), ymax=c(Inf,Inf), fill="Gleason \u2264 6"),
                alpha=0.25) +
    geom_segment(aes(x=prev_biopsies, xend=prev_biopsies,
                     y=-rep(Inf,length(prev_biopsies)), 
                     yend=rep(Inf,length(prev_biopsies)),
                     linetype=rep("Previous Biopsies", length(prev_biopsies)))) + 
    geom_vline(xintercept = current_visit_time, linetype='dashed') +
    geom_line(data=schedule_df,
              aes(x=biopsy_times, y=schedule_id, group=schedule_id), color=THEME_COLOR,
              linetype='dotted') +
    geom_label(data=schedule_df,
               aes(x=biopsy_times, y=schedule_id, label='B', group=schedule_id),
               fill=THEME_COLOR,
               color='white',
               size=7) +
    geom_label(aes(x=c(0, current_visit_time), 
                   y=rep(tail(ybreaks,1)+1,2),
                   label=c(startAS_label, current_visit_label)), 
               fill=WARNING_COLOR, color='white', size=6, label.r = unit(0.25,units = "lines")) +
    scale_fill_manual(name="", values = SUCCESS_COLOR) + 
    scale_linetype_manual(name="", values = c("solid"))+
    theme_bw() + xlab("Follow-up time (years)") +
    scale_x_continuous(breaks = xTicks, labels=xlabels,
                       limits = c(0, MAX_FOLLOW_UP))+
    scale_y_continuous(breaks = ybreaks, labels = ylabels, 
                       limits = c(ybreaks[1]-0.5, tail(ybreaks,1) + 1.5)) +
    ylab("Biopsy Schedule") + 
    theme(text = element_text(size = FONT_SIZE),
          axis.text.y = element_text(size = FONT_SIZE),
          axis.text.x = element_text(size = FONT_SIZE,
                                     angle = 30, hjust = 1),
          legend.text = element_text(size=FONT_SIZE),
          legend.position = "bottom", legend.direction = "horizontal")
  
  if(!is.null(latest_biopsy_label)){
    pp = pp + geom_label(aes(x=latest_survival_time, y=sum(range(ybreaks))/2, 
                             label=latest_biopsy_label), 
                         fill=WARNING_COLOR, color='white', size=6, label.r = unit(0.25,units = "lines"))
  }
  
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

biopsyDelayGaugeGraph = function(delay, schedule, max_delay=DELAY_GAUGE_MAX){
  prop_delay = delay/max_delay
  
  gauge_color = colorRamp(c(SUCCESS_COLOR, WARNING_COLOR, DANGER_COLOR))(prop_delay)
  gauge_color = rgb(gauge_color[1], gauge_color[2], gauge_color[3], maxColorValue = 255)
  
  gauge_ticks_colors = sapply(seq(0,1,0.25), FUN = function(prop){
    col=colorRamp(c(SUCCESS_COLOR, WARNING_COLOR, DANGER_COLOR))(prop)
    rgb(col[1], col[2], col[3], maxColorValue = 255)
  })
  
  risk_label = paste0("\n\nDelay: ", round(delay)," months\nSchedule: ", schedule)
  
  riskGauge = ggplot(data = NULL, 
                     aes(ymax = prop_delay, ymin = 0, xmax = 2, xmin = 1, 
                         fill="Risk")) +
    geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="white", color=gauge_color) +
    geom_rect() +
    geom_segment(aes(x=2.0, xend=2.1, y=0, yend=0), color=gauge_ticks_colors[1])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.25, yend=0.25), color=gauge_ticks_colors[2])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.5, yend=0.5), color=gauge_ticks_colors[3])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.75, yend=0.75), color=gauge_ticks_colors[4])+
    geom_segment(aes(x=2.0, xend=2.1, y=1, yend=1), color=gauge_ticks_colors[5])+
    geom_text(aes(x = 2.25, y = 0, label = 0), size=4, color=gauge_ticks_colors[1]) +
    geom_text(aes(x = 2.25, y = 0.25, label = round(max_delay * 0.25,1)), size=4, color=gauge_ticks_colors[2]) +
    geom_text(aes(x = 2.25, y = 0.5, label = round(max_delay * 0.5,1)), size=4, color=gauge_ticks_colors[3]) +
    geom_text(aes(x = 2.25, y = 0.75, label = round(max_delay * 0.75,1)), size=4, color=gauge_ticks_colors[4]) +
    geom_text(aes(x = 2.25, y = 1, label = max_delay), size=4, color=gauge_ticks_colors[5]) +
    geom_text(aes(x = 0, y = 0, label = risk_label), 
              color=gauge_color, size=6) +
    coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2.5)) + ylim(c(0,2)) +
    geom_segment(aes(x=0, xend=1, y=prop_delay, yend=prop_delay), 
                 color=gauge_color,
                 arrow = arrow(length = unit(0.2,"cm")))+
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

riskGaugeGraph = function(mean_risk_prob, date="", danger_color_threshold = 0.25){
  
  if(mean_risk_prob > danger_color_threshold){
    gauge_color = DANGER_COLOR
  }else{
    gauge_color = colorRamp(c(SUCCESS_COLOR, WARNING_COLOR, DANGER_COLOR))(mean_risk_prob/danger_color_threshold)
    gauge_color = rgb(gauge_color[1], gauge_color[2], gauge_color[3], maxColorValue = 255)
  }
  
  risk_label = paste0("\n\n\nCumulative-risk: ", round(mean_risk_prob*100),"%\n on ",date)
  
  gauge_ticks_colors = sapply(seq(0,1,0.25), FUN = function(prop){
    if(prop > danger_color_threshold){
      return(DANGER_COLOR)
    }else{
      col=colorRamp(c(SUCCESS_COLOR, WARNING_COLOR, DANGER_COLOR))(prop/danger_color_threshold)
      return(rgb(col[1], col[2], col[3], maxColorValue = 255))
    }
  })
  
  riskGauge = ggplot(data = NULL, 
                     aes(ymax = mean_risk_prob, ymin = 0, xmax = 2, xmin = 1, 
                         fill="Risk")) +
    geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="white", color=gauge_color) +
    geom_rect() +
    geom_segment(aes(x=2.0, xend=2.1, y=0, yend=0), color=gauge_ticks_colors[1])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.25, yend=0.25), color=gauge_ticks_colors[2])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.5, yend=0.5), color=gauge_ticks_colors[3])+
    geom_segment(aes(x=2.0, xend=2.1, y=0.75, yend=0.75), color=gauge_ticks_colors[4])+
    geom_segment(aes(x=2.0, xend=2.1, y=1, yend=1), color=gauge_ticks_colors[5])+
    geom_text(aes(x = 2.35, y = 0, label = "0%"), size=5, color=gauge_ticks_colors[1]) +
    geom_text(aes(x = 2.35, y = 0.25, label = "25%"), size=5, color=gauge_ticks_colors[2]) +
    geom_text(aes(x = 2.35, y = 0.5, label = "50%"), size=5, color=gauge_ticks_colors[3]) +
    geom_text(aes(x = 2.35, y = 0.75, label = "75%"), size=5, color=gauge_ticks_colors[4]) +
    geom_text(aes(x = 2.35, y = 1, label = "100%"), size=5, color=gauge_ticks_colors[5]) +
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


