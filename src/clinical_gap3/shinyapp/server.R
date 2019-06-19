library(shiny)
library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)
library(plyr)

shinyServer(function(input, output, session) {
  
  #making a reactive value to track changes in patient data
  patientCounter <- reactiveVal(0)
  biopsyCounter <- reactiveVal(0)
  
  resetPatientCache = function(){
    patient_cache <<- list(P_ID = -1, data=NULL)
  }
  
  recalculateBiopsySchedules = function(){
    patient_data = patient_cache$patient_data
    biopsy_schedules = compareSchedules(patient_data, patient_cache$current_visit_time,
                                        patient_cache$latest_survival_time,
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
    biopsyCounter(biopsyCounter() + 1)
  }
  
  setPatientDataInCache = function(patient_data){
    patient_data = patient_data[patient_data$year_visit<=MAX_FOLLOW_UP,]
    
    if(max(patient_data$gleason_sum, na.rm=T) > 6){
      patient_data$gleason_sum[patient_data$gleason_sum %in% c(7,8,9,10)] = NA
    }
    
    patient_cache <<- list(P_ID = patient_data$P_ID[1],
                           patient_data = patient_data,
                           dom_diagnosis = patient_data$dom_diagnosis[1],
                           latest_survival_time = max(patient_data$year_visit[!is.na(patient_data$gleason_sum)]),
                           current_visit_time = max(patient_data$year_visit))
    
    set.seed(patient_cache$P_ID)
    latest_survival_time = max(patient_data$year_visit[!is.na(patient_data$gleason_sum)])
    
    patient_cache$SURV_CACHE_TIMES <<- c(patient_cache$latest_survival_time, 
                                         seq(patient_cache$latest_survival_time+1/365, MAX_FOLLOW_UP, length.out = CACHE_SIZE-1))
    patient_cache$PSA_CACHE_TIMES <<- c(seq(0, patient_cache$latest_survival_time - 1/365, length.out = 20), patient_cache$SURV_CACHE_TIMES)
    
    pred_res = getExpectedFutureOutcomes(mvJoint_psa_time_scaled, patient_data, 
                                         latest_survival_time, 
                                         survival_predict_times = patient_cache$SURV_CACHE_TIMES[-1],
                                         psa_predict_times = patient_cache$PSA_CACHE_TIMES, 
                                         psaDist = "Tdist", M = M, addRandomError = T)
    
    patient_cache$PSA_CACHE_FULL <<- pred_res$predicted_psa
    patient_cache$SURV_CACHE_MEAN <<- c(1,rowMeans(pred_res$predicted_surv_prob))
    patient_cache$SURV_CACHE_FULL <<- rbind(rep(1, M), pred_res$predicted_surv_prob)
    
    #slider on tab 2
    sliderLabel = paste0("Time (years) since latest visit in ", 
                         getHumanReadableDate(patient_cache$dom_diagnosis + patient_cache$current_visit_time*YEAR_DIVIDER, abbreviated = T))
    
    max_slider = round_any(MAX_FOLLOW_UP - patient_cache$current_visit_time, accuracy=0.5)
    updateSliderInput(session, "risk_pred_time", value = 0,
                      max = max_slider, label = sliderLabel)
    
    recalculateBiopsySchedules()
    
    patientCounter(patientCounter() + 1)
  }
  
  observeEvent(input$load_pat1, {
    setPatientDataInCache(demo_pat_list[[1]])
  })
  
  observeEvent(input$load_pat2, {
    setPatientDataInCache(demo_pat_list[[2]])
  })
  
  observeEvent(input$load_pat3, {
    setPatientDataInCache(demo_pat_list[[3]])
  })
  
  observeEvent(input$load_pat4, {
    setPatientDataInCache(demo_pat_list[[4]])
  })
  
  observeEvent(input$patientFile,{
    inFile <- input$patientFile
    
    if (!is.null(inFile)){
      patient_data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                              sep = input$sep, quote = input$quote)
      setPatientDataInCache(patient_data)
    }else{
      resetPatientCache()
    }
  })
  
  observeEvent(input$year_gap_biopsy,{
    if(patientCounter()>0 & is.numeric(input$year_gap_biopsy)){
      recalculateBiopsySchedules()
    }
  })
  
  #For the first tab panel's first table
  output$table_obs_data <- renderTable({
    if(patientCounter()>0){
      patient_data = patient_cache$patient_data
      
      first_visit_date = getHumanReadableDate(patient_cache$dom_diagnosis)
      last_biopsy_date = getHumanReadableDate(patient_cache$dom_diagnosis + patient_cache$latest_survival_time*YEAR_DIVIDER)
      latest_visit_date = getHumanReadableDate(patient_cache$dom_diagnosis + patient_cache$current_visit_time*YEAR_DIVIDER)
      nr_biopsies = sum(!is.na(patient_data$gleason_sum))
      nr_visits = nrow(patient_data)
      age = round(patient_data$age[1])
      latest_gleason = tail(patient_data$gleason_sum[!is.na(patient_data$gleason_sum)],1)
      
      return(data.frame("Data"=c("ID", "Age (years)", 
                                 "First Visit", "Latest Visit", "Total Visits", 
                                 "Latest Biopsy", "Latest Gleason","Total Biopsies"),
                        "Value"=c(patient_cache$P_ID, age,
                                  first_visit_date,latest_visit_date, nr_visits,
                                  last_biopsy_date, latest_gleason,nr_biopsies)))
      
    }else{
      return(NULL)
    }
  })
  
  output$graph_obs_psa <- renderPlot({
    if(patientCounter()>0){
      sampled_psa_indices = c(1:20, seq(22, length(patient_cache$PSA_CACHE_TIMES), by = 7))
      return(psaObsDataGraph(patient_cache$patient_data,
                             patient_cache$PSA_CACHE_TIMES[sampled_psa_indices],
                             patient_cache$PSA_CACHE_FULL[sampled_psa_indices,]))
    }else{
      return(NULL)
    }
  })
  
  output$cum_risk_gauge = renderPlot({
    if(patientCounter()>0){
      futureTime = patient_cache$current_visit_time + input$risk_pred_time
      riskProb = 1 - patient_cache$SURV_CACHE_MEAN[which.min(abs(patient_cache$SURV_CACHE_TIMES - futureTime))]
      date = getHumanReadableDate(patient_cache$dom_diagnosis + futureTime*YEAR_DIVIDER)
      
      return(riskGaugeGraph(riskProb, date))
    }else{
      return(NULL)
    }
  })
  
  output$biopsy_schedule_graph = renderPlot({
    if(patientCounter()>0 & biopsyCounter()>0 & !is.null(input$selected_schedules)){
      return(biopsyScheduleGraph(patient_cache$biopsy_schedule_plotDf,
                                 as.numeric(input$selected_schedules),
                                 patient_cache$dom_diagnosis,
                                 patient_cache$current_visit_time,
                                 patient_cache$latest_survival_time))
    }else{
      return(NULL)
    }
  }) 
  
  output$biopsy_total_delay_graph = renderPlot({
    if(patientCounter()>0 & biopsyCounter()>0 & !is.null(input$selected_schedules)){
      biopsy_total_graph = biopsyTotalGraph(patient_cache$biopsy_total_delay_plotDf,
                                            as.numeric(input$selected_schedules))
      biopsy_delay_graph = biopsyDelayGraph(patient_cache$biopsy_total_delay_plotDf,
                                            as.numeric(input$selected_schedules))
      return(ggarrange(biopsy_total_graph, biopsy_delay_graph,
                       align = "h", widths = c(1.15,1)))
    }else{
      return(NULL)
    }
  })
  
  resetPatientCache()
})
