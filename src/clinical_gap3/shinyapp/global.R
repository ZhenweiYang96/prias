load("mvJoint_psa_time_scaled_light.Rdata")
load("demo_pat_list.Rdata")
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
SCHEDULES = c("5% Risk", "10% Risk", "15% Risk",
              "PRIAS", "Yearly", "Every 2 Years")
DELAY_GAUGE_MAX = 24

DANGER_COLOR = "#c71c22"
THEME_COLOR = "#2fa4e7"
THEME_COLOR_DARK = "#033a6f"
SUCCESS_COLOR = "#73A839"
WARNING_COLOR="orange"

getHumanReadableDate = function(spss_date, abbreviated=F){
  format(as.POSIXct(spss_date, origin = SPSS_ORIGIN_DATE), 
         format = ifelse(abbreviated, DATE_PRINT_FORMAT_ABBREVIATED, DATE_PRINT_FORMAT))
}

psaObsDataGraph = function(data, dom_diagnosis, current_visit_time, 
                           latest_survival_time,
                           pred_times, pred_log2psaplus1){
  biopsy_times = data$year_visit[!is.na(data$gleason_sum)]
  
  pred_mean_psa = rowMeans(2^pred_log2psaplus1 - 1)
  
  first_visit_date = getHumanReadableDate(dom_diagnosis, abbreviated = T)
  last_biopsy_date = getHumanReadableDate(dom_diagnosis + latest_survival_time*YEAR_DIVIDER, abbreviated = T)
  current_visit_date = getHumanReadableDate(dom_diagnosis + current_visit_time*YEAR_DIVIDER, abbreviated = T)
  
  x_points = seq(0, current_visit_time, length.out = 5)
  x_labels = as.character(round(x_points,1))
  
  latest_biopsy_point_done = F
  x_labels[1] = "0\n(Start AS"
  if(x_points[1]==latest_survival_time){
    x_labels[1] = paste0(x_labels[1], " &\nLatest biopsy\n", first_visit_date, ")")
    latest_biopsy_point_done = T
  }else{
    x_labels[1] = paste0(x_labels[1],"\n",first_visit_date,")")
  }
  
  x_labels[5] = paste(round(current_visit_time,1),"\n(Current visit")
  if(x_points[5]==latest_survival_time){
    x_labels[5] = paste0(x_labels[5], " &\nLatest biopsy\n", current_visit_date, ")")
    latest_biopsy_point_done = T
  }else{
    x_labels[5] = paste0(x_labels[5],"\n",current_visit_date,")")
  }
  
  if(latest_biopsy_point_done==F){
    pp = c(2:4)[which.min(abs(x_points[c(2,3,4)] - latest_survival_time))]
    x_points[[pp]] = latest_survival_time
    if((current_visit_time - latest_survival_time) < 1 | latest_survival_time<1){
      x_labels[[pp]] = paste(round(latest_survival_time,1),"\n(Latest      \nbiopsy    \n", last_biopsy_date, ")    ")
    }else{
      x_labels[[pp]] = paste(round(latest_survival_time,1),"\n(Latest biopsy\n", last_biopsy_date, ")")
    }
    
    #readjust other points
    if(pp==3){
      x_points[2] = latest_survival_time/2
      x_labels[2] = round(x_points[2],1)
      
      x_points[4] = (latest_survival_time + current_visit_time)/2
      x_labels[4] = round(x_points[4],1)
    }else if(pp==4){
      x_points[2:3] = seq(0, latest_survival_time, length.out = 4)[2:3]
      x_labels[2:2] = round(x_points[2:3],1)
    }else if(pp==2){
      x_points[3:4] = seq(latest_survival_time, current_visit_time, length.out = 4)[2:3]
      x_labels[3:4] = round(x_points[3:4],1)
    }
  }
  
  plot = ggplot() +
    geom_ribbon(aes(x=c(0, latest_survival_time),
                    ymin=-c(Inf,Inf), ymax=c(Inf,Inf), fill="Gleason ≤ 6"),
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
    scale_fill_manual(name="", values = SUCCESS_COLOR) + 
    scale_shape_manual(name="", values = 16)+
    scale_linetype_manual(name="", values = c("solid","dotted"))+
    scale_x_continuous(breaks = x_points, labels = x_labels,
                       limits = c(0, 1.1*current_visit_time))+
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
                               latest_survival_time,
                               prev_biopsies){
  
  schedule_df = schedule_df[schedule_df$schedule_id %in% selected_schedule_ids,]
  
  if(current_visit_time == latest_survival_time){
    xTicks = c(0, seq(current_visit_time, MAX_FOLLOW_UP, length.out = 5))  
  }else{
    xTicks = c(0, latest_survival_time, seq(current_visit_time, MAX_FOLLOW_UP, length.out = 4))
  }
  
  xlabels = as.character(round(xTicks,1))
  xTicks_spps_dates = sapply(start_dom + xTicks*YEAR_DIVIDER, getHumanReadableDate, abbreviated=T)
  xlabels[1] = paste0("0\n(Start AS\n", xTicks_spps_dates[1], ")")
  
  if(current_visit_time == latest_survival_time){
    xlabels[2] = paste0(xlabels[2],"\n(Current visit\n& Latest biopsy\n",xTicks_spps_dates[2], ")")
    xlabels[3:6] = sapply(3:6, function(i){paste0(xlabels[i],"\n(",xTicks_spps_dates[i], ")")})
  }else{
    if((current_visit_time - latest_survival_time) < 1 | latest_survival_time < 1){
      xlabels[2] = paste0(xlabels[2],"\n(Latest      \n biopsy      \n",xTicks_spps_dates[2], ")      ")
    }else{
      xlabels[2] = paste0(xlabels[2],"\n(Latest biopsy\n",xTicks_spps_dates[2], ")")
    }
    
    xlabels[3] = paste0(xlabels[3],"\n(Current visit\n",xTicks_spps_dates[3], ")")
    xlabels[4:6] = sapply(4:6, function(i){paste0(xlabels[i],"\n(",xTicks_spps_dates[i], ")")})
  }
  
  ybreaks = unique(schedule_df$schedule_id)
  ylabels = unique(schedule_df$Schedule)
  
  pp = ggplot() +
    geom_ribbon(aes(x=c(0, latest_survival_time),
                    ymin=-c(Inf,Inf), ymax=c(Inf,Inf), fill="Gleason ≤ 6"),
                alpha=0.25) +
    geom_segment(aes(x=prev_biopsies, xend=prev_biopsies,
                     y=-rep(Inf,length(prev_biopsies)), 
                     yend=rep(Inf,length(prev_biopsies)),
                     linetype=rep("Previous Biopsies", length(prev_biopsies)))) + 
    geom_vline(xintercept = current_visit_time, linetype='dashed') +
    geom_label(data=schedule_df,
               aes(x=biopsy_times, y=schedule_id, label='B', group=schedule_id), color=THEME_COLOR,
               size=7) +
    geom_line(data=schedule_df,
              aes(x=biopsy_times, y=schedule_id, group=schedule_id), color=THEME_COLOR,
              linetype='dotted') +
    scale_fill_manual(name="", values = SUCCESS_COLOR) + 
    scale_linetype_manual(name="", values = c("solid"))+
    theme_bw() + xlab("Follow-up time (years)") +
    scale_x_continuous(breaks = xTicks, labels=xlabels,
                       limits = c(0, MAX_FOLLOW_UP))+
    scale_y_continuous(breaks = ybreaks, labels = ylabels, 
                       limits = c(ybreaks[1]-0.5, tail(ybreaks,1) + 0.5)) +
    ylab("Biopsy Schedule") + 
    theme(text = element_text(size = FONT_SIZE + 2),
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

riskGaugeGraph = function(mean_risk_prob, date=""){
  
  gauge_color = colorRamp(c(SUCCESS_COLOR, WARNING_COLOR, DANGER_COLOR))(mean_risk_prob)
  gauge_color = rgb(gauge_color[1], gauge_color[2], gauge_color[3], maxColorValue = 255)
  
  risk_label = paste0("\n\n\nCumulative-risk: ", round(mean_risk_prob*100,2),"%\n on ",date)
  
  gauge_ticks_colors = sapply(seq(0,1,0.25), FUN = function(prop){
    col=colorRamp(c(SUCCESS_COLOR, WARNING_COLOR, DANGER_COLOR))(prop)
    rgb(col[1], col[2], col[3], maxColorValue = 255)
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


