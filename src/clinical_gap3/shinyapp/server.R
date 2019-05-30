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

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
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
    
    pid = data$P_ID[1]
    first_visit_date = getHumanReadableDate(data$dom_diagnosis[1])
    last_biopsy_date = getHumanReadableDate(max(data$dom_visit[!is.na(data$gleason_sum)]))
    latest_visit_date = getHumanReadableDate(max(data$dom_visit))
    nr_biopsies = sum(!is.na(data$gleason_sum))
    nr_visits = nrow(data)
    age = round(data$age[1])
    latest_gleason = tail(data$gleason_sum[!is.na(data$gleason_sum)],1)
  
    data_to_show = data.frame("Data"=c("ID", "Age (years)", 
                                               "First Visit", "Last Visit", "Total Visits", 
                                               "Last Biopsy", "Last Gleason","Total Biopsies"),
                              "Value"=c(pid, age,
                                        first_visit_date,latest_visit_date, nr_visits,
                                        last_biopsy_date, latest_gleason,nr_biopsies))
    
    return(data_to_show)
  })
  
  output$graph_obs_psa <- renderPlot({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$patientFile
    
    if (is.null(inFile)){
      return(NULL)
    }
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data = data[data$year_visit<=MAX_FOLLOW_UP,]
   
    return(psaObsDataGraph(data))
  })
  
  output$cum_risk_gauge = renderPlot({
    inFile <- input$patientFile
    
    if (is.null(inFile)){
      return(NULL)
    }
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data = data[data$year_visit<=MAX_FOLLOW_UP,]
    
    futureVisitTime = data$year_max_followup[1] + input$risk_pred_time
    return(riskGaugeGraph(data, futureVisitTime))
  })
  
  output$table_biopsy_options = renderTable({
    inFile <- input$patientFile
    
    if (is.null(inFile))
      return(NULL)
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data = data[data$year_visit<=MAX_FOLLOW_UP,]
    
    horizon = input$horizon_biopsy
    min_gap_biopsy = input$year_gap_biopsy
    
    output$biopsy_options_graph = renderPlot(getBiopsyOptionsGraph(options_df))
    
    return(NULL)
  })
})
