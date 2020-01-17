library(shiny)
library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)
library(plyr)
library(xlsx)

shinyServer(function(input, output, session) {
  
  options(warn = -1)
  #####################
  #Functions related to patient cache
  #####################
  #making a reactive value to track changes in patient data
  patientCounter <- reactiveVal(0)
  biopsyCounter <- reactiveVal(0)
  
  first_time <<- TRUE
  
  resetPatientCache = function(){
    patient_cache <<- list(P_ID = -1, patient_data=NULL)
  }
  
  recalculateBiopsySchedules = function(){
    pat_data = patient_cache$patient_data
    cur_visit_time = patient_cache$current_visit_time
    latest_survival_time = patient_cache$latest_survival_time
    min_biosy_gap = input$month_gap_biopsy/12
    
    set.seed(2019)
    schedule_5perc = getRiskBasedSchedule(object = mvJoint_psa_time_scaled, patient_data = pat_data,
                                          cur_visit_time = cur_visit_time, 
                                          latest_survival_time = latest_survival_time,
                                          risk_threshold = 0.05, min_biopsy_gap = min_biosy_gap, M = M,
                                          horizon = MAX_FOLLOW_UP)
    schedule_10perc = getRiskBasedSchedule(object = mvJoint_psa_time_scaled, patient_data = pat_data,
                                           cur_visit_time = cur_visit_time, 
                                           latest_survival_time = latest_survival_time,
                                           risk_threshold = 0.1, min_biopsy_gap = min_biosy_gap, M = M,
                                           horizon = MAX_FOLLOW_UP)
    schedule_15perc = getRiskBasedSchedule(object = mvJoint_psa_time_scaled, patient_data = pat_data,
                                           cur_visit_time = cur_visit_time, 
                                           latest_survival_time = latest_survival_time,
                                           risk_threshold = 0.15, min_biopsy_gap = min_biosy_gap, M = M,
                                           horizon = MAX_FOLLOW_UP)
    
    schedule_annual = getFixedSchedule(cur_visit_time = cur_visit_time, 
                                       latest_survival_time = latest_survival_time,
                                       min_biopsy_gap = min_biosy_gap, biopsy_frequency = 1, 
                                       horizon = MAX_FOLLOW_UP)
    schedule_biennial = getFixedSchedule(cur_visit_time = cur_visit_time, 
                                         latest_survival_time = latest_survival_time,
                                         min_biopsy_gap = min_biosy_gap, biopsy_frequency = 2, 
                                         horizon = MAX_FOLLOW_UP)
    
    
    schedule_prias = getPRIASSchedule(object = mvJoint_psa_time_scaled, patient_data = pat_data,
                                      cur_visit_time = cur_visit_time, latest_survival_time = latest_survival_time,
                                      min_biopsy_gap = 1, M = M, horizon = MAX_FOLLOW_UP)
    
    consequences_list = lapply(list(schedule_5perc, schedule_10perc,schedule_15perc,
                                    schedule_annual,schedule_biennial, schedule_prias),
                               FUN = function(x){
                                 set.seed(2019);
                                 getConsequences(object=mvJoint_psa_time_scaled, 
                                                 patient_data=pat_data,
                                                 cur_visit_time = cur_visit_time, 
                                                 proposed_biopsy_times = x,
                                                 latest_survival_time = latest_survival_time, M=M,
                                                 horizon=MAX_FOLLOW_UP)
                               })
    
    expected_delays = sapply(consequences_list, "[[", "expected_delay")*12
    practical_biopsy_times = lapply(consequences_list, "[[", "practical_biopsy_times")
    total_biopsies = sapply(practical_biopsy_times, length)
    
    patient_cache$biopsy_schedule_plotDf <<- data.frame(schedule_id = rep(1:length(SCHEDULES), total_biopsies),
                                                        Schedule=rep(SCHEDULES, total_biopsies), 
                                                        biopsy_times=unlist(practical_biopsy_times))  
    
    patient_cache$biopsy_total_delay_plotDf <<- data.frame(schedule_id = 1:length(SCHEDULES),
                                                           schedule=SCHEDULES,
                                                           expected_delay=expected_delays,
                                                           total_biopsies = total_biopsies, check.names = F)
    biopsyCounter(biopsyCounter() + 1)
  }
  
  setPatientDataInCache = function(patient_data){
    showModal(pleaseWaitModal)
    
    start_date = as.Date(patient_data$start_date, tryFormats = c("%d-%m-%Y"))
    visit_date = as.Date(patient_data$visit_date, tryFormats = c("%d-%m-%Y"))
    
    patient_data$log2psaplus1 = log(patient_data$psa + 1, base = 2)
    patient_data$dom_diagnosis = as.numeric(difftime(start_date, SPSS_ORIGIN_DATE, units='secs'))
    patient_data$year_visit = (as.numeric(difftime(visit_date, SPSS_ORIGIN_DATE, units='secs')) - patient_data$dom_diagnosis)/YEAR_DIVIDER
    
    patient_data = patient_data[patient_data$year_visit<=MAX_FOLLOW_UP,]
    if(max(patient_data$gleason_sum, na.rm=T) > 6){
      patient_data$gleason_sum[patient_data$gleason_sum %in% c(7,8,9,10)] = NA
    }
    
    patient_cache <<- list(P_ID = patient_data$P_ID[1],
                           patient_data = patient_data,
                           dom_diagnosis = patient_data$dom_diagnosis[1],
                           latest_survival_time = max(patient_data$year_visit[!is.na(patient_data$gleason_sum)]),
                           current_visit_time = max(patient_data$year_visit))
    
    cur_visit_time_in_secs = patient_cache$dom_diagnosis + patient_cache$current_visit_time*YEAR_DIVIDER
    max_visit_time_in_secs = patient_cache$dom_diagnosis + MAX_FOLLOW_UP * YEAR_DIVIDER
    
    set.seed(2019)
    latest_survival_time = max(patient_data$year_visit[!is.na(patient_data$gleason_sum)])
    
    patient_cache$SURV_CACHE_TIMES <<- unique(c(seq(patient_cache$current_visit_time, MAX_FOLLOW_UP, STEP_CUMRISK_SLIDER),MAX_FOLLOW_UP))
    patient_cache$PSA_CACHE_TIMES <<- seq(0, MAX_FOLLOW_UP, length.out = 50)
    
    pred_res = getExpectedFutureOutcomes(mvJoint_psa_time_scaled, patient_data, 
                                         latest_survival_time, 
                                         survival_predict_times = patient_cache$SURV_CACHE_TIMES,
                                         psa_predict_times = patient_cache$PSA_CACHE_TIMES, 
                                         psaDist = "Tdist", M = M, addRandomError = F)
    
    patient_cache$PSA_CACHE_FULL <<- pred_res$predicted_psa
    if(latest_survival_time == MAX_FOLLOW_UP){
      patient_cache$SURV_CACHE_MEAN = 1
    }else{
      patient_cache$SURV_CACHE_MEAN <<- rowMeans(pred_res$predicted_surv_prob, na.rm = T)
      if(patient_cache$current_visit_time == latest_survival_time){
        patient_cache$SURV_CACHE_MEAN <<- c(1, patient_cache$SURV_CACHE_MEAN)
      }
    }
    
    #slider on tab 2
    updateSliderInput(session, "risk_pred_time", value = as.POSIXct(cur_visit_time_in_secs, origin = SPSS_ORIGIN_DATE),
                      max = as.POSIXct(max_visit_time_in_secs, origin = SPSS_ORIGIN_DATE), 
                      min = as.POSIXct(cur_visit_time_in_secs, origin = SPSS_ORIGIN_DATE),
                      step = NULL,
                      label = "Choose a time to predict cumulative-risk of reclassification:",
                      timeFormat = "%b %Y")
    
    recalculateBiopsySchedules()
    
    patientCounter(patientCounter() + 1)
    removeModal()
  }
  
  ######################
  # All modals
  ######################
  
  #modal to show processing
  pleaseWaitModal = modalDialog(title = "Please wait", "Analyzing patient data.", 
                                footer = NULL,size = "s", fade = F)
  
  #modal to show cohort change
  cohort_change_modal_PRIAS = modalDialog(
    title="Selected Cohort: PRIAS",
    tags$p("In this cohort, risk predictions and personalized biopsy schedules can only be made for the first 6 years of follow-up. Patient data beyond this limit is automatically trimmed."),
    tags$p("Removing any previously loaded patient data."),
    footer = modalButton("Ok"),size = "m", fade = F, easyClose = T
  )
  
  cohort_change_modal_Toronto = modalDialog(
    title="Selected Cohort: Toronto AS",
    tags$p("In this cohort, risk predictions and personalized biopsy schedules can only be made for the first 8 years of follow-up. Patient data beyond this limit is automatically trimmed."),
    tags$p("Removing any previously loaded patient data."),
    footer = modalButton("Ok"),size = "m", fade = F, easyClose = T
  )
  
  cohort_change_modal_Hopkins = modalDialog(
    title="Selected Cohort: Johns Hopkins AS",
    tags$p("In this cohort, risk predictions and personalized biopsy schedules can only be made for the first 7 years of follow-up. Patient data beyond this limit is automatically trimmed."),
    tags$p("Removing any previously loaded patient data."),
    footer = modalButton("Ok"),size = "m", fade = F, easyClose = T
  )
  
  cohort_change_modal_MSKCC = modalDialog(
    title="Selected Cohort: MSKCC AS",
    tags$p("In this cohort, risk predictions and personalized biopsy schedules can only be made for the first 6 years of follow-up. Patient data beyond this limit is automatically trimmed."),
    tags$p("Removing any previously loaded patient data."),
    footer = modalButton("Ok"),size = "m", fade = F, easyClose = T
  )
  
  cohort_change_modal_KCL = modalDialog(
    title="Selected Cohort: KCL (London) AS",
    tags$p("In this cohort, risk predictions and personalized biopsy schedules can only be made for the first 3 years of follow-up. Patient data beyond this limit is automatically trimmed."),
    tags$p("Removing any previously loaded patient data."),
    footer = modalButton("Ok"),size = "m", fade = F, easyClose = T
  )
  
  cohort_change_modal_MUSIC = modalDialog(
    title="Selected Cohort: MUSIC AS",
    tags$p("In this cohort, risk predictions and personalized biopsy schedules can only be made for the first 2 years of follow-up. Patient data beyond this limit is automatically trimmed."),
    tags$p("Removing any previously loaded patient data."),
    footer = modalButton("Ok"),size = "m", fade = F, easyClose = T
  )
  
  cohort_change_modal_UCSF = modalDialog(
    title="Selected Cohort: UCSF AS",
    tags$p("In this cohort, risk predictions and personalized biopsy schedules can only be made for the first 8.5 years of follow-up. Patient data beyond this limit is automatically trimmed."),
    tags$p("Removing any previously loaded patient data."),
    footer = modalButton("Ok"),size = "m", fade = F, easyClose = T
  )
  
  MODAL_MAPPING = list("PRIAS"=cohort_change_modal_PRIAS,
                       "Toronto"=cohort_change_modal_Toronto,
                       "Hopkins"=cohort_change_modal_Hopkins,
                       "MSKCC"=cohort_change_modal_MSKCC,
                       "MUSIC"=cohort_change_modal_MUSIC,
                       "KCL"=cohort_change_modal_KCL,
                       "UCSF"=cohort_change_modal_UCSF)
  
  #modal for manual entry
  manual_entry_modal = modalDialog(
    tags$h3("Patient Data Manual Entry Form"),
    tags$hr(),
    tags$h4(tags$span("Dates are in dd-mm-yyyy format.", class='label label-warning')),
    numericInput("manual_age", label = "Enter patient age (years)", value = 60, width = "80%"),
    textInput("manual_dom_diagnosis", label="Enter date (dd-mm-yyyy format) of low-grade prostate cancer diagnosis:",
              value = "01-01-2019", width = "80%"),
    textInput("manual_biopsy_times", 
              label = "Enter dates (dd-mm-yyyy format) of previous biopsies with Gleason grade 1 (Gleason 3+3). Separate dates by comma.",
              value = "01-01-2019, 02-01-2020, 03-06-2021",
              width = "80%"),
    textInput("manual_psa_times", 
              label = "Enter dates (dd-mm-yyyy format) of previous PSA measurements. Separate dates by comma.",
              value = "01-01-2019, 04-06-2019, 05-01-2020, 08-06-2020, 02-01-2021, 03-06-2021, 02-01-2022, 09-06-2022",
              width = "80%"),
    textInput("manual_psa_values", 
              label = "Enter PSA values (ng/mL). Separate them by comma.",
              value = "5.7, 3.2, 12, 8.5, 15, 21.7, 25, 20.3", width = "80%"),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("ok_manual_entry", "OK", class='btn-primary')
    )
  )
  
  #Modal to show example data
  exampleDataModal = modalDialog(
    tags$h3("Example Patient Excel Data"),
    tags$hr(),
    tags$h4(tags$span("All column names are case sensitive.", class='label label-danger')),
    tags$h4(tags$span("Missing values should be left blank.", class='label label-warning')),
    tags$h4(tags$span("Dates are in dd-mm-yyyy format.", class='label label-warning')),
    tableOutput('example_data_in_modal'),
    tags$hr(),
    tags$span("Description", class="lead"),
    tags$br(),
    tags$span("age", class='label label-primary'),
    tags$span(" is the age (years) of the patient when patient started AS. Missing values are not allowed."),
    tags$br(),
    tags$span("start_date", class='label label-primary'),
    tags$span(" is the date (dd-mm-yyyy format) on which patient started AS. Missing values are not allowed."),
    tags$br(),
    tags$span("visit_date", class='label label-primary'),
    tags$span(" is the date (dd-mm-yyyy format) of follow-up on which either PSA was measured or a biopsy was conducted. Missing values are not allowed. It should never be before the start date."),
    tags$br(),
    tags$span("psa", class='label label-primary'),
    tags$span(" is the PSA (ng/mL) at the follow-up time. Missing values should be left blank."),
    tags$br(),
    tags$span("gleason_grade", class='label label-primary'),
    tags$span(" is the ISUP Gleason grade group (lowest 1, highest 5) at the follow-up time. Missing values should be left blank."),
    tags$br(),
    tags$hr(),
    downloadButton("download_example_data2", "Download Example File", class='btn-success'),
    footer = tagList(modalButton("OK"))
  )
  
  output$example_data_in_modal = renderTable(EXAMPLE_DF)
  exampleDataDownloadHandler = downloadHandler(
    filename = "example_dataset.xlsx",
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    content = function(file) {
      write.xlsx(x=EXAMPLE_DF, file = file,row.names = FALSE, col.names = T,
                 append = F, showNA = F)
    }
  )
  output$download_example_data = exampleDataDownloadHandler
  output$download_example_data2 = exampleDataDownloadHandler
  
  ################################
  # Events for buttons on sidebar panel
  ################################
  observeEvent(input$cohort, {
    mvJoint_psa_time_scaled <<- models[[input$cohort]]
    MAX_FOLLOW_UP <<- MAX_FOLLOW_UP_MAPPING[input$cohort]
    CURRENT_COHORT_NAME <<- names(COHORT_MAPPING[which(COHORT_MAPPING==input$cohort)])
    
    if(!first_time){
      showModal(MODAL_MAPPING[[input$cohort]])
    } else{
      first_time <<- FALSE
    }
    
    resetPatientCache()
    patientCounter(patientCounter() + 1)
  })
  
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
  
  observeEvent(input$load_manual,{
    showModal(manual_entry_modal)
  })
  
  observeEvent(input$ok_manual_entry, {
    manual_biopsy_date = as.Date(unlist(strsplit(input$manual_biopsy_times, ",")),tryFormats = c("%d-%m-%Y"))
    manual_psa_date = as.Date(unlist(strsplit(input$manual_psa_times, ",")), tryFormats = c("%d-%m-%Y"))
    manual_psa_values = as.numeric(unlist(strsplit(input$manual_psa_values, ",")))
    
    visit_date = sort(unique(c(manual_biopsy_date, manual_psa_date)))
    
    psa = rep(NA, length(visit_date))
    gleason_sum = rep(NA, length(visit_date))
    psa[visit_date %in% manual_psa_date] = manual_psa_values
    gleason_sum[visit_date %in% manual_biopsy_date] = 6
    
    manual_pat_data = data.frame('P_ID'=-10, age = input$manual_age,
                                 start_date=input$manual_dom_diagnosis,
                                 visit_date=format(visit_date, "%d-%m-%Y"),
                                 psa = psa, gleason_sum=gleason_sum)
    
    removeModal()
    
    setPatientDataInCache(manual_pat_data)
  }) 
  
  observeEvent(input$show_example_data, {
    showModal(exampleDataModal)
  })
  
  observeEvent(input$patientFile,{
    inFile <- input$patientFile
    
    if (!is.null(inFile)){
      patient_data = read.xlsx(inFile$datapath,sheetIndex = 1,
                               as.data.frame = T, header = T)
      
      patient_data$P_ID = -15
      patient_data$psa[patient_data$psa %in% "NA"] = NA
      patient_data$gleason_grade[patient_data$gleason_grade %in% "NA"] = NA
      patient_data$psa = as.numeric(as.character(patient_data$psa))
      patient_data$gleason_grade = as.numeric(as.character(patient_data$gleason_grade))
      patient_data$gleason_sum[!is.na(patient_data$gleason_grade) & 
                                 patient_data$gleason_grade==1] = 6
      patient_data$gleason_sum[!is.na(patient_data$gleason_grade) & 
                                 patient_data$gleason_grade > 1] = 7
      
      setPatientDataInCache(patient_data)
    }else{
      resetPatientCache()
    }
  })
  
  #############################
  # Output on Tab 1 showing patient data
  #############################
  output$table_obs_data <- renderTable({
    if(patientCounter()>0 & !is.null(patient_cache$patient_data)){
      patient_data = patient_cache$patient_data
      
      first_visit_date = getHumanReadableDate(patient_cache$dom_diagnosis)
      last_biopsy_date = getHumanReadableDate(patient_cache$dom_diagnosis + patient_cache$latest_survival_time*YEAR_DIVIDER)
      current_visit_date = getHumanReadableDate(patient_cache$dom_diagnosis + patient_cache$current_visit_time*YEAR_DIVIDER)
      nr_biopsies = sum(!is.na(patient_data$gleason_sum))
      nr_visits = nrow(patient_data)
      age = round(patient_data$age[1])
      latest_gleason = tail(patient_data$gleason_sum[!is.na(patient_data$gleason_sum)],1)
      
      if(latest_gleason <= 6){
        latest_gleason_grade = "1 (Gleason 3+3)"
      } 
      
      return(data.frame("Data"=c("Age (years)", 
                                 "First Visit", "Current Visit", "Total Visits", 
                                 "Latest Biopsy", "Latest Gleason Grade","Total Biopsies"),
                        "Value"=c(age,
                                  first_visit_date,current_visit_date, nr_visits,
                                  last_biopsy_date, latest_gleason_grade, nr_biopsies)))
      
    }else{
      return(NULL)
    }
  })
  
  output$graph_obs_psa <- renderPlot({
    if(patientCounter()>0 & !is.null(patient_cache$patient_data)){
      return(psaObsDataGraph(patient_cache$patient_data,
                             patient_cache$dom_diagnosis,
                             patient_cache$current_visit_time,
                             patient_cache$latest_survival_time,
                             patient_cache$PSA_CACHE_TIMES,
                             patient_cache$PSA_CACHE_FULL))
    }else{
      return(NULL)
    }
  })
  
  #############################
  # Output on Tab 2 showing patient risk
  #############################
  output$table_obs_data_cumrisk <- renderTable({
    if(patientCounter()>0 & !is.null(patient_cache$patient_data)){
      patient_data = patient_cache$patient_data
      
      first_visit_date = getHumanReadableDate(patient_cache$dom_diagnosis)
      last_biopsy_date = getHumanReadableDate(patient_cache$dom_diagnosis + patient_cache$latest_survival_time*YEAR_DIVIDER)
      current_visit_date = getHumanReadableDate(patient_cache$dom_diagnosis + patient_cache$current_visit_time*YEAR_DIVIDER)
      nr_biopsies = sum(!is.na(patient_data$gleason_sum))
      nr_visits = nrow(patient_data)
      age = round(patient_data$age[1])
      latest_gleason = tail(patient_data$gleason_sum[!is.na(patient_data$gleason_sum)],1)
      
      if(latest_gleason <= 6){
        latest_gleason_grade = "1 (Gleason 3+3)"
      } 
      
      return(data.frame("Data"=c("Age (years)", 
                                 "First Visit", "Current Visit", "Total Visits", 
                                 "Latest Biopsy", "Latest Gleason Grade", "Total Biopsies"),
                        "Value"=c(age,
                                  first_visit_date,current_visit_date, nr_visits,
                                  last_biopsy_date, latest_gleason_grade, nr_biopsies)))
      
    }else{
      return(NULL)
    }
  })
  
  output$cum_risk_gauge = renderPlot({
    if(patientCounter()>0 & input$risk_pred_time!=0 & !is.null(patient_cache$patient_data)){
      futureTime = as.numeric((difftime(input$risk_pred_time, SPSS_ORIGIN_DATE, units='secs') - patient_cache$dom_diagnosis)/YEAR_DIVIDER)
      
      riskProb = 1 - patient_cache$SURV_CACHE_MEAN[which.min(abs(patient_cache$SURV_CACHE_TIMES - futureTime))]
      date = getHumanReadableDate(patient_cache$dom_diagnosis + futureTime*YEAR_DIVIDER)
      
      return(riskGaugeGraph(riskProb, date))
    }else{
      return(NULL)
    }
  })
  
  #############################
  # Output on Tab 3 showing biopsy recommendation
  #############################
  observeEvent(input$month_gap_biopsy,{
    if(patientCounter()>0 & is.numeric(input$month_gap_biopsy) & !is.null(patient_cache$patient_data)){
      showModal(pleaseWaitModal)
      recalculateBiopsySchedules()
      removeModal()
    }
  })
  
  decisionRenderer = function(decision_label_num){
    renderUI({
      if(patientCounter()>0 & biopsyCounter()>0 & !is.null(patient_cache$patient_data) &
         !is.null(input$selected_schedules) & length(input$selected_schedules)> (decision_label_num-1)){
        
        schedule_id = as.numeric(input$selected_schedules[decision_label_num])
        
        schedule_name = paste0(names(SCHEDULES_MAPPING)[schedule_id], ":")
        first_biopsy_time = patient_cache$biopsy_schedule_plotDf$biopsy_times[patient_cache$biopsy_schedule_plotDf$schedule_id==schedule_id][1]
        
        first_biopsy_date = getHumanReadableDate(patient_cache$dom_diagnosis + first_biopsy_time*YEAR_DIVIDER)
        #class = paste("label", ifelse(first_biopsy_time > patient_cache$current_visit_time, "label-primary", "label-success"))
        class = "label label-default"
        
        return(tags$div(tags$span(schedule_name, class=class), 
                        tags$span(first_biopsy_date, class='text-default'),
                        class='lead'))
      }else{
        return(NULL)
      }
    })
  }
  output$decision1 = decisionRenderer(1)
  output$decision2 = decisionRenderer(2)
  output$decision3 = decisionRenderer(3)
  output$decision4 = decisionRenderer(4)
  output$decision5 = decisionRenderer(5)
  output$decision6 = decisionRenderer(6)
  
  output$biopsy_schedule_graph = renderPlot({
    if(patientCounter()>0 & biopsyCounter()>0 & !is.null(patient_cache$patient_data) & !is.null(input$selected_schedules)){
      return(biopsyScheduleGraph(patient_cache$biopsy_schedule_plotDf,
                                 as.numeric(input$selected_schedules),
                                 patient_cache$dom_diagnosis,
                                 patient_cache$current_visit_time,
                                 patient_cache$latest_survival_time,
                                 patient_cache$patient_data$year_visit[!is.na(patient_cache$patient_data$gleason_sum)]))
    }else{
      return(NULL)
    }
  }) 
  
  output$biopsy_delay_gauge_graph = renderPlot({
    if(patientCounter()>0 & biopsyCounter()>0 & !is.null(patient_cache$patient_data) & !is.null(input$selected_schedules)){
      max_delay = max(DELAY_GAUGE_MAX, 
                      ceiling(patient_cache$biopsy_total_delay_plotDf$expected_delay[as.numeric(input$selected_schedules)]))
      
      plotDf = patient_cache$biopsy_total_delay_plotDf
      
      delayGaugeList = lapply(1:length(input$selected_schedules), FUN = function(x){
        roww = plotDf$schedule_id == as.numeric(input$selected_schedules[x])
        biopsyDelayGaugeGraph(plotDf$expected_delay[roww], plotDf$schedule[roww], max_delay = max_delay)
      })
      
      return(ggarrange(plotlist = delayGaugeList, nrow=2, ncol=3))
    }else{
      return(NULL)
    }
  })
  
  output$delay_explanation_plot = renderPlot({
    ggplot() + geom_ribbon(aes(x=c(3.25, 4.5), ymin=-Inf, ymax=Inf), fill=DANGER_COLOR, alpha=0.25) +
      geom_label(aes(x=3.875, y=0.8, label='15 months time\ndelay in detecting\nGleason reclassification'),size=5)+
      geom_segment(aes(x=c(0,1,2), y=c(-Inf, -Inf, -Inf), xend=c(0,1,2), yend=c(0.5, 0.5, 0.5)),
                   color=c(WARNING_COLOR, rep(SUCCESS_COLOR, 2)))+
      geom_vline(xintercept = 4.5, color=DANGER_COLOR) + 
      geom_vline(xintercept = 3.25, linetype='dashed', color=DANGER_COLOR) + 
      geom_label(aes(x=3.25, y=0.5, label = "True time of\nGleason\nreclassification"),
                 size=6, color='white',
                 fill=c(DANGER_COLOR)) +
      geom_label(aes(x=c(0,1,2,4.5), y=rep(0.5,4),
                     label = c("Start AS", rep("Biopsy\nGleason grade 1",2), "Biopsy Gleason\nreclassification\ndetected")),
                 size=6, color=c('white', rep(SUCCESS_COLOR,2), DANGER_COLOR),
                 fill=c(WARNING_COLOR, rep('white',3))) +
      xlab("Time of biopsy visits") + 
      scale_x_continuous(breaks = c(0,1,2,3.25,4.5), 
                         labels = c("Jan 2005","Jan 2006", "Jan 2007",
                                    "Apr 2008", "Jul 2009"), limits=c(-0.1,5)) + 
      ylim(0,1) + 
      theme_bw() +
      theme(text = element_text(size=FONT_SIZE + 4), 
            axis.text.x = element_text(size=FONT_SIZE + 4), 
            axis.title.x= element_text(size=FONT_SIZE + 4), 
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), axis.ticks.y = element_blank())
  })
  
  #Finally before we start the server, we reset the patient cache
  resetPatientCache()
})
