#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shinythemes)
library(shiny)
library(shinyalert)
library(ggplot2)
library(plotly)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = shinytheme("cerulean"),
                  useShinyalert(),
                  
                  tags$head(
                    tags$link(rel = "stylesheet", type = "text/css", href = "main.css")
                  ),
                  
                  # Application title
                  #tags$div(class="text-center bg-dark",
                  tags$h1(style="font-size: 40px; margin: 0 0 0 0; padding: 25px 0px 25px 0px;",
                          class="text-center bg-primary no-margin","Biopsy Recommender for Prostate Cancer Patients in Active Surveillance"),
                  
                  tags$hr(),
                  
                  sidebarLayout(
                    # Sidebar with a slider input for number of bins 
                    sidebarPanel(
                      tags$label("Load demo patient data"),
                      tags$br(), 
                      actionButton("load_pat1", "Patient 1", class="btn-default"),
                      tags$br(), tags$br(),
                      actionButton("load_pat2", "Patient 2", class="btn-default"),
                      tags$br(), tags$br(),
                      actionButton("load_pat3", "Patient 3", class="btn-default"),
                      tags$br(), tags$br(),
                      actionButton("load_pat4", "Patient 4", class="btn-default"),
                      tags$hr(),
                      tags$label("Enter patient data manually"),
                      actionButton("load_manual", "Manual Entry", class="btn-primary"),
                      tags$hr(),
                      fileInput('patientFile', 'Load patient data (Excel file)',
                                accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                      actionButton("show_example_data", "Show Example File", class='btn-info'),
                      tags$br(),tags$br(),
                      downloadButton("download_example_data", "Download Example File", class='btn-success'),
                      width=2
                    ),
                    
                    # Show a plot of the generated distribution
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel("Patient Data", value="patient_data", tags$br(),
                                           fluidRow(column(5, tags$div(class="panel panel-default",
                                                                       tags$div(class="panel-heading",
                                                                                "Patient History"),
                                                                       tableOutput("table_obs_data"))),
                                                    column(7, plotOutput("graph_obs_psa")))
                                  ),
                                  tabPanel("Risk of Gleason ≥ 7", tags$br(),
                                           sliderInput("risk_pred_time", "Time (years) since latest visit:",
                                                       0, 10, 0, step = 0.5, animate = animationOptions(interval = 2000, 
                                                                                                        loop = FALSE)),
                                           plotOutput("cum_risk_gauge", width = "400px"),
                                           textOutput("gauge_date")),
                                  tabPanel("Biopsy recommendation", 
                                           tags$br(), 
                                           tags$span("Schedule Description", class='lead'),
                                           tags$br(),
                                           tags$span("Fixed: Yearly", class='label label-default'),
                                           tags$span("schedules a biopsy every year."),
                                           tags$br(),
                                           tags$span("Fixed: Every 2 years", class='label label-default'),
                                           tags$span("schedules a biopsy every 2 years."),
                                           tags$br(),
                                           tags$span("Fixed: PRIAS", class='label label-default'),
                                           tags$span("schedules a biopsy at year 1, 4, 7 and 10. It switches to yearly biopsies if PSA doubling time is between 0 and 10."),
                                           tags$br(),
                                           tags$span("Personalized: 5% risk", class='label label-success'),
                                           tags$span("schedules a biopsy at a visit if a patient's cumulative risk of Gleason ≥ 7 is above 5%."),
                                           tags$br(),
                                           tags$span("Personalized: 10% risk", class='label label-success'),
                                           tags$span("schedules a biopsy at a visit if a patient's cumulative risk of Gleason ≥ 7 is above 10%."),
                                           tags$br(),
                                           tags$span("Personalized: 15% risk", class='label label-success'),
                                           tags$span("schedules a biopsy at a visit if a patient's cumulative risk of Gleason ≥ 7 is above 15%."),
                                           tags$br(),
                                           tags$hr(),
                                           fluidRow(column(6, numericInput("year_gap_biopsy",
                                                                           "Minimum gap (years) between biopsies",
                                                                           step=0.5, min=0, value=1)),
                                                    column(6, selectizeInput('selected_schedules',
                                                                             label = "Select Biopsy Schedule", 
                                                                             choices = c("Personalized: 5% risk"=1,
                                                                                         "Personalized: 10% risk"=2,
                                                                                         "Personalized: 15% risk"=3,
                                                                                         "Fixed: PRIAS"=4,
                                                                                         "Fixed: Yearly"=5,
                                                                                         "Fixed: Every 2 years"=6),
                                                                             selected = c(2,4,5),
                                                                             multiple=TRUE))),
                                           tags$span("Recommended Date of Next Biopsy", class='lead'),
                                           fluidRow(column(4, htmlOutput("decision1")),
                                                    column(4, htmlOutput("decision2")),
                                                    column(4, htmlOutput("decision3"))),
                                           fluidRow(column(4, htmlOutput("decision4")),
                                                    column(4, htmlOutput("decision5")),
                                                    column(4, htmlOutput("decision6"))),
                                           tags$br(), tags$br(),
                                           tags$span("Recommended Dates of Future Biopsies", class='lead'),
                                           plotOutput("biopsy_schedule_graph"),
                                           tags$br(), tags$br(),
                                           tags$span("Expected Time Delay (months) to Detect Gleason ≥ 7", class='lead'),
                                           plotOutput("biopsy_delay_gauge_graph"))
                                           
                                  )
                      )
                    )
                  )
)
