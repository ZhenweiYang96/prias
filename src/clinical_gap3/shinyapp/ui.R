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
  
  # Application title
  #tags$div(class="text-center bg-dark",
  tags$h1(style="font-size: 40px; margin: 0 0 0 0; padding: 25px 0px 25px 0px;",
          class="text-center bg-primary no-margin","Biopsy Recommender for Prostate Cancer Patients in Active Surveillance"),
  
  tags$hr(),

  sidebarLayout(
    # Sidebar with a slider input for number of bins 
    sidebarPanel(
      
      #fileInput('RDfile', 'Load the R Workspace with the fitted joint model',
      #          accept = NULL),
      
      fileInput('patientFile', 'Load patient data (CSV file)',
                accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      radioButtons('dec', 'Decimal', c(Dot = '.', Comma = ','), 
                   selected='.'),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"'),
      width=2
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Patient Data", value="patient_data", tags$br(), 
                           splitLayout(cellWidths = c("18%","5%", "65%"),
                             tags$div(class="panel panel-default",
                                                tags$div(class="panel-heading",
                                                         "Patient History"),
                                                tableOutput("table_obs_data")),
                             tags$div(), plotOutput("graph_obs_psa", height="500px"))
                           ),
                  tabPanel("Biopsy recommendation", tags$br(), 
                           splitLayout(verticalLayout(
                             titlePanel("Cumulative-risk of Gleason Upgrade"),
                             plotOutput("plot_cum_risk")),
                    #sliderInput("risk_pred_time", "Years from now:", 
                    #            0, 5, 0, step = 0.5),
                    #plotOutput("cum_risk_gauge")),
                          
                    verticalLayout(titlePanel("Biopsy Options"),
                                   numericInput("year_gap_biopsy",
                                                "Minimum gap (years) between biopsies",
                                                step=0.5, min=0, value=1),
                                   tags$div(class="panel panel-default",
                                            tags$div(class="panel-heading",
                                                     "Biopsy Options"),
                                            tableOutput("table_biopsy_options")),
                                   plotOutput("biopsy_options_graph"))
                  ))
      )
    )
  ))
)
  