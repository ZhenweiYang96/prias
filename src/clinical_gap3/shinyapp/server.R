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
  
  resetPatientCache = function(){
    patient_cache <<- list(P_ID = -1)
  }
  
  setSurvDataInCache = function(pat_data){
    patient_cache$P_ID <<- pat_data$P_ID[1]
    if(is.null(patient_cache$SURV_CACHE_FULL)){
      set.seed(patient_cache$P_ID)
      latest_survival_time = max(pat_data$year_visit[!is.na(pat_data$gleason_sum)])
      
      patient_cache$SURV_CACHE_TIMES <<- c(latest_survival_time, 
                                           seq(latest_survival_time+1/365, MAX_FOLLOW_UP, length.out = CACHE_SIZE-1))
      patient_cache$PSA_CACHE_TIMES <<- c(seq(0, latest_survival_time - 1/365, length.out = 20), patient_cache$SURV_CACHE_TIMES)
      pred_res = getExpectedFutureOutcomes(mvJoint_psa_time_scaled, pat_data, 
                                           latest_survival_time, 
                                           survival_predict_times = patient_cache$SURV_CACHE_TIMES[-1],
                                           psa_predict_times = patient_cache$PSA_CACHE_TIMES, 
                                           psaDist = "Tdist", M = M, addRandomError = T)
      patient_cache$PSA_CACHE_FULL <<- pred_res$predicted_psa
      patient_cache$SURV_CACHE_MEAN <<- c(1,rowMeans(pred_res$predicted_surv_prob))
      patient_cache$SURV_CACHE_FULL <<- rbind(rep(1, M), pred_res$predicted_surv_prob)
      rm(pred_res)
    }
  }
  
  output$table_obs_data <- renderTable({
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
        resetPatientCache()
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
        setSurvDataInCache(data)
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
      resetPatientCache()
      
      sliderLabel = paste0("Time (years) since latest visit in ", 
                           getHumanReadableDate(max(data$dom_visit), abbreviated = T))
      max_slider = round_any(MAX_FOLLOW_UP - max(data$year_visit), accuracy=0.5)
      updateSliderInput(session, "risk_pred_time", value = 0,
                        max = max_slider, label = sliderLabel)
    }
    
    patient_cache$P_ID <<- data$P_ID[1]
    if(is.null(patient_cache$SURV_CACHE_FULL)){
      setSurvDataInCache(data)
    }
    
    futureTime = max(data$year_visit) + input$risk_pred_time
    riskProb = 1 - patient_cache$SURV_CACHE_MEAN[which.min(abs(patient_cache$SURV_CACHE_TIMES - futureTime))]
    date = getHumanReadableDate(data$dom_diagnosis[1] + futureTime*YEAR_DIVIDER)
    
    return(riskGaugeGraph(riskProb, date))
  })
  
  output$biopsy_schedule_graph = renderPlot({
    inFile <- input$patientFile
    
    if (is.null(inFile)){
      return(NULL)
    }
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data = data[data$year_visit<=MAX_FOLLOW_UP,]
    
    if(patient_cache$P_ID!=data$P_ID[1]){
      resetPatientCache()
    }
    
    patient_cache$P_ID <<- data$P_ID[1]
    if(is.null(patient_cache$SURV_CACHE_FULL)){
      setSurvDataInCache(data)
    }
    
    if(is.null(patient_cache$biopsy_schedule_plotDf)){
      biopsy_schedules = compareSchedules(data, max(data$year_visit),
                                          max(data$year_visit[!is.na(data$gleason_sum)]),
                                          risk_thresholds = c(0.05, 0.1, 0.15),
                                          input$year_gap_biopsy,
                                          patient_cache$SURV_CACHE_TIMES, patient_cache$PSA_CACHE_TIMES,
                                          patient_cache$SURV_CACHE_FULL, patient_cache$PSA_CACHE_FULL)
      
      schedules = c("5% Risk", "10% Risk", "15% Risk",
                    "PRIAS", "Yearly", "Every 2 Years")
      expected_delays = sapply(biopsy_schedules$schedules, "[[", "expected_delay") * 12
      total_biopsies = sapply(biopsy_schedules$schedules, "[[", "total_biopsies")
      
      biopsy_times = do.call('c', lapply(biopsy_schedules$schedules, "[[", "biopsy_times"))
      patient_cache$biopsy_schedule_plotDf <<- data.frame(schedule_id = rep(1:length(schedules), total_biopsies),
                                                          Schedule=rep(schedules, total_biopsies),
                                                          biopsy_times)
      
      patient_cache$biopsy_total_delay_plotDf <<- data.frame(schedule_id = 1:length(schedules),
                                                             schedule=schedules,
                                                             expected_delay=expected_delays,
                                                             total_biopsies = total_biopsies, check.names = F)
    }
    
    if(is.null(input$selected_schedules)){
      output$biopsy_total_delay_graph = renderPlot(NULL)
      return(NULL)
    }else{
      biopsy_schedule_graph = biopsyScheduleGraph(patient_cache$biopsy_schedule_plotDf,
                                                  as.numeric(input$selected_schedules),
                                                  data$dom_diagnosis[1],
                                                  max(data$year_visit),
                                                  max(data$year_visit[!is.na(data$gleason_sum)]))
      
      biopsy_total_graph = biopsyTotalGraph(patient_cache$biopsy_total_delay_plotDf,
                                                              as.numeric(input$selected_schedules))
      biopsy_delay_graph = biopsyDelayGraph(patient_cache$biopsy_total_delay_plotDf,
                                                              as.numeric(input$selected_schedules))
      output$biopsy_total_delay_graph = renderPlot(ggarrange(biopsy_total_graph, biopsy_delay_graph, 
                                                  align = "h", widths = c(1.15,1)))
      return(biopsy_schedule_graph)
    }
  })
  
  resetPatientCache()
})
