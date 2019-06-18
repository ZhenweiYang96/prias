#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#]
library(shiny)
library(JMbayes)
library(splines)
library(ggplot2)
library(plotly)
library(DT)
library(ggpubr)
library(plyr)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  patient_cache <<- list(P_ID = -1)
  
  output$table_obs_data <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$patientFile
    
    if (is.null(inFile))
      return(NULL)
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data = data[data$year_visit<=MAX_FOLLOW_UP,]
    
    #patient already had an event
    if(max(data$gleason_sum, na.rm = T)>6){
      return(NULL)
    }
    
    if(is.null(patient_cache$data_to_show) | patient_cache$P_ID!=data$P_ID[1]){
      if(patient_cache$P_ID!=data$P_ID[1]){
        patient_cache <<- list(P_ID = -1)
      }
      
      first_visit_date = getHumanReadableDate(data$dom_diagnosis[1])
      last_biopsy_date = getHumanReadableDate(max(data$dom_visit[!is.na(data$gleason_sum)]))
      latest_visit_date = getHumanReadableDate(max(data$dom_visit))
      nr_biopsies = sum(!is.na(data$gleason_sum))
      nr_visits = nrow(data)
      age = round(data$age[1])
      latest_gleason = tail(data$gleason_sum[!is.na(data$gleason_sum)],1)
      
      patient_cache$P_ID <<- data$P_ID[1]
      
      if(is.null(patient_cache$SURV_CACHE_FULL)){
        set.seed(patient_cache$P_ID)
        latest_survival_time = max(data$year_visit[!is.na(data$gleason_sum)])
        
        patient_cache$SURV_CACHE_TIMES <<- c(latest_survival_time, 
                                             seq(latest_survival_time+1/365, MAX_FOLLOW_UP, length.out = CACHE_SIZE-1))
        patient_cache$PSA_CACHE_TIMES <<- c(seq(0, latest_survival_time - 1/365, length.out = 20), patient_cache$SURV_CACHE_TIMES)
        pred_res = getExpectedFutureOutcomes(mvJoint_psa_time_scaled, data, 
                                             latest_survival_time, 
                                             survival_predict_times = patient_cache$SURV_CACHE_TIMES[-1],
                                             psa_predict_times = patient_cache$PSA_CACHE_TIMES, 
                                             psaDist = "Tdist", M = M, addRandomError = T)
        patient_cache$PSA_CACHE_FULL <<- pred_res$predicted_psa
        patient_cache$SURV_CACHE_MEAN <<- c(1,rowMeans(pred_res$predicted_surv_prob))
        patient_cache$SURV_CACHE_FULL <<- rbind(rep(1, M), pred_res$predicted_surv_prob)
        rm(pred_res)
      }
      
      patient_cache$data_to_show <<- data.frame("Data"=c("ID", "Age (years)", 
                                                         "First Visit", "Latest Visit", "Total Visits", 
                                                         "Latest Biopsy", "Latest Gleason","Total Biopsies"),
                                                "Value"=c(patient_cache$P_ID, age,
                                                          first_visit_date,latest_visit_date, nr_visits,
                                                          last_biopsy_date, latest_gleason,nr_biopsies))
      
      sampled_psa_indices = c(1:20, seq(22, length(patient_cache$PSA_CACHE_TIMES), by = 7))
      patient_cache$psa_obs_data_graph <<- psaObsDataGraph(data, 
                                                           patient_cache$PSA_CACHE_TIMES[sampled_psa_indices], 
                                                           patient_cache$PSA_CACHE_FULL[sampled_psa_indices,])
    }
    
    output$graph_obs_psa = renderPlot(patient_cache$psa_obs_data_graph)
    
    return(patient_cache$data_to_show)
  })
  
  output$cum_risk_gauge = renderPlot({
    inFile <- input$patientFile
    
    if (is.null(inFile)){
      return(NULL)
    }
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data = data[data$year_visit<=MAX_FOLLOW_UP,]
    
    if(patient_cache$P_ID!=data$P_ID[1]){
      patient_cache <<- list(P_ID = -1)
      sliderLabel = paste0("Time (years) since latest visit in ", 
                     getHumanReadableDate(max(data$dom_visit), abbreviated = T))
      max_slider = round_any(MAX_FOLLOW_UP - max(data$year_visit), accuracy=0.5)
      updateSliderInput(session, "risk_pred_time", value = 0,
                        max = max_slider, label = sliderLabel)
    }
    
    patient_cache$P_ID <<- data$P_ID[1]
    if(is.null(patient_cache$SURV_CACHE_FULL)){
      set.seed(patient_cache$P_ID)
      latest_survival_time = max(data$year_visit[!is.na(data$gleason_sum)])
      
      patient_cache$SURV_CACHE_TIMES <<- c(latest_survival_time, 
                                           seq(latest_survival_time+1/365, MAX_FOLLOW_UP, length.out = CACHE_SIZE-1))
      patient_cache$PSA_CACHE_TIMES <<- c(seq(0, latest_survival_time - 1/365, length.out = 20), patient_cache$SURV_CACHE_TIMES)
      pred_res = getExpectedFutureOutcomes(mvJoint_psa_time_scaled, data, 
                                           latest_survival_time, 
                                           survival_predict_times = patient_cache$SURV_CACHE_TIMES[-1],
                                           psa_predict_times = patient_cache$PSA_CACHE_TIMES, 
                                           psaDist = "Tdist", M = M, addRandomError = T)
      patient_cache$PSA_CACHE_FULL <<- pred_res$predicted_psa
      patient_cache$SURV_CACHE_MEAN <<- c(1,rowMeans(pred_res$predicted_surv_prob))
      patient_cache$SURV_CACHE_FULL <<- rbind(rep(1, M), pred_res$predicted_surv_prob)
      rm(pred_res)
    }
    
    futureTime = max(data$year_visit) + input$risk_pred_time
    riskProb = 1 - patient_cache$SURV_CACHE_MEAN[which.min(abs(patient_cache$SURV_CACHE_TIMES - futureTime))]
    date = getHumanReadableDate(data$dom_diagnosis[1] + futureTime*YEAR_DIVIDER)
    
    return(riskGaugeGraph(riskProb, date))
  })
  
  output$biopsy_options_graph = renderPlot({
    # biopsy_schedules = compareSchedules(data, max(data$year_visit),
    #                                     max(data$year_visit[!is.na(data$gleason_sum)]),
    #                                     risk_thresholds = c(0.05, 0.1, 0.15),
    #                                     input$year_gap_biopsy, 
    #                                     patient_cache$SURV_CACHE_TIMES, patient_cache$PSA_CACHE_TIMES,
    #                                     patient_cache$SURV_CACHE_FULL, patient_cache$PSA_CACHE_FULL)
    # 
    # schedules = c("5% Risk", "10% Risk", "15% Risk",
    #               "PRIAS", "Yearly", "Every 2 Years")
    # expected_delays = sapply(biopsy_schedules$schedules, "[[", "expected_delay") * 12
    # total_biopsies = sapply(biopsy_schedules$schedules, "[[", "total_biopsies")
    # 
    # biopsy_times = do.call('c', lapply(biopsy_schedules$schedules, "[[", "biopsy_times"))
    # plotDf = data.frame(Schedule=rep(schedules, total_biopsies),
    #                     biopsy_times)
    # 
    # patient_cache$schedule_tab <<- data.frame("Schedule"=schedules,
    #                                           "Delay (months)"=expected_delays,
    #                                           "Total Biopsies" = total_biopsies, check.names = F)
    # 
    # xTicks = seq(max(data$year_visit),
    #              MAX_FOLLOW_UP, length.out = 4)
    # xTicks_spps_dates = data$dom_diagnosis[1] + xTicks*YEAR_DIVIDER
    # xlabels = sapply(xTicks_spps_dates, getHumanReadableDate, abbreviated=T)
    # xlabels[1] = paste0(getHumanReadableDate(xTicks_spps_dates[1]),"\n(Current Visit)")
    # 
    # patient_cache$biopsy_times_plot <<- ggplot(plotDf) +
    #   geom_label(aes(x=biopsy_times, y=Schedule, label='B'), color=THEME_COLOR,
    #              size=7) +
    #   geom_line(aes(x=biopsy_times, y=Schedule), color=THEME_COLOR,
    #             linetype='dotted') +
    #   theme_bw() + xlab("Proposed Dates of Biopsy") +
    #   scale_x_continuous(breaks = xTicks, labels=xlabels,
    #                      limits = c(data$year_max_followup[1], MAX_FOLLOW_UP))+
    #   ylab("Biopsy Schedule") + theme(axis.text = element_text(size = FONT_SIZE),
    #                                   axis.title = element_text(size = FONT_SIZE),
    #                                   legend.position = "bottom", legend.direction = "horizontal")
    # 
    # output$biopsy_options_graph = renderPlot(patient_cache$biopsy_times_plot)
    # output$table_biopsy_options = renderTable(patient_cache$schedule_tab)
  })
})
