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
                  
                  tags$br(),
                  
                  sidebarLayout(
                    # Sidebar with a slider input for number of bins 
                    sidebarPanel(
                      selectInput("cohort", "Choose AS cohort:", 
                                  choices=COHORT_MAPPING,
                                  selected = "PRIAS"),
                      hr(),
                      tags$h4("Load Patient Data"),
                      tags$br(), 
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
                      fileInput('patientFile', 'or Load patient file (Excel format)',
                                accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                      actionButton("show_example_data", "Show Example File", class='btn-info'),
                      tags$br(),tags$br(),
                      downloadButton("download_example_data", "Download Example File", class='btn-success'),
                      #tags$hr(),
                      #tags$label("or Enter patient data manually"),
                      #actionButton("load_manual", "Manual Entry", class="btn-primary"),
                      width=2
                    ),
                    
                    # Show a plot of the generated distribution
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel("Patient Summary", value="patient_data", tags$br(),
                                           fluidRow(column(5, tags$div(class="panel panel-default",
                                                                       tags$div(class="panel-heading",
                                                                                "Observed Data"),
                                                                       tableOutput("table_obs_data")))),
                                           fluidRow(column(7, tags$span("Observed and model fitted PSA", class='lead'),
                                                           tags$br(),
                                                           plotOutput("graph_obs_psa")),
                                                    column(5, tags$span("Future cumulative-risk of Gleason upgrading", class='lead'),
                                                           tags$br(),
                                                           sliderInput("risk_pred_time", 
                                                                       "Choose a time to predict cumulative-risk of Upgrading:",
                                                                       min = as.POSIXct(0, origin=SPSS_ORIGIN_DATE), 
                                                                       max = as.POSIXct(MAX_FOLLOW_UP*YEAR_DIVIDER, origin=SPSS_ORIGIN_DATE), 
                                                                       value = as.POSIXct(0, origin=SPSS_ORIGIN_DATE), step = YEAR_DIVIDER*STEP_CUMRISK_SLIDER, timeFormat="%b %Y",
                                                                       animate = animationOptions(interval = 2000, loop = FALSE), width="400px"),
                                                           plotOutput("cum_risk_gauge", width = "400px")))
                                  ),
                                  
                                  tabPanel("Biopsy Recommendation", 
                                           tags$br(),
                                           tags$div(style="display: none", numericInput("month_gap_biopsy",
                                                                                        "Time gap (months) between consecutive biopsies",
                                                                                        step=6, min=0, value=12)),
                                           fluidRow(column(4, tags$div(tags$h4("Personalized Schedule"), tags$hr(),
                                                                       tags$strong("Description:"),
                                                                       tags$span("schedules a biopsy at a visit if a patient's cumulative risk of Gleason upgrading is above a risk threshold. This risk threshold can be chosen using the slider below. A gap of 12 months is maintained between consecutive biopsies."),
                                                                       tags$br(),
                                                                       tags$strong("Default threshold:"),
                                                                       tags$span("By default a patient-specific threshold is chosen automatically as described in Documentation point number 5."),
                                                                       tags$br(),tags$br(),
                                                                       sliderInput("risk_threshold", 
                                                                                   "Choose a risk threshold (%):",
                                                                                   min = 0, 
                                                                                   max = 100, value=10,
                                                                                   step = 1),
                                                                       tags$strong("Timing of planned biopsies"),
                                                                       tags$div(class="panel panel-default",
                                                                                tableOutput("table_planned_personalized")),
                                                                       tags$strong("Expected time delay (months) in detecting Gleason upgrading:"),
                                                                       plotOutput("biopsy_delay_gauge_graph_personalized"), 
                                                                       class="roundedCorner")),
                                                    column(4, tags$div(tags$h4("Yearly Schedule"), tags$hr(),
                                                                       tags$strong("Description:"),
                                                                       tags$span("schedules a biopsy every year."),
                                                                       tags$br(),tags$br(),
                                                                       tags$strong("Timing of planned biopsies"),
                                                                       tags$div(class="panel panel-default",
                                                                                tableOutput("table_planned_yearly")),
                                                                       tags$strong("Expected time delay (months) in detecting Gleason upgrading"),
                                                                       plotOutput("biopsy_delay_gauge_graph_yearly"),
                                                                       class="roundedCorner")),
                                                    column(4, tags$div(tags$h4("PRIAS Schedule"), tags$hr(),
                                                                       tags$strong("Description:"),
                                                                       tags$span("schedules a biopsy at year 1, 4, 7 and 10. It switches to yearly biopsies if PSA doubling time is between 0 and 10."),
                                                                       tags$br(), tags$br(),
                                                                       tags$strong("Timing of planned biopsies"),
                                                                       tags$div(class="panel panel-default",
                                                                                tableOutput("table_planned_prias")),
                                                                       tags$strong("Expected time delay (months) in detecting Gleason upgrading"),
                                                                       plotOutput("biopsy_delay_gauge_graph_prias"),
                                                                       class="roundedCorner")))
                                  ),
                                  tabPanel("Documentation", 
                                           tags$h3(tags$span("Upgrading is defined as increase in Gleason grade group from 1 (Gleason 3+3) to 2 (Gleason 3+4) or higher.", class='label label-info')),
                                           tags$h3(tags$span("The risk predictions, biopsy schedules, and expected delays in detection of Gleason Upgrading are only indicative.", class='label label-info')),
                                           tags$h3(tags$span("Clinical consultation is recommended.", class='label label-warning')),
                                           tags$hr(),
                                           tags$h2(tags$span("1. Calculation of cumulative-risk of Gleason Upgrading", class='label label-success')),
                                           tags$br(),
                                           tags$img(src="jmpred.png", width="70%"),
                                           tags$hr(),
                                           tags$h2(tags$span("2. Prediction joint model", class='label label-success')),
                                           tags$span("We fitted a ", class='lead'),
                                           tags$a("joint model for time-to-event and longitudinal data [1]",target="_blank",
                                                  href="https://www.jstatsoft.org/index.php/jss/article/view/v072i07/v72i07.pdf", class='lead'),
                                           tags$span(" using the ", class='lead'),
                                           tags$a("R package JMbayes", target="_blank", href='https://cran.r-project.org/web/packages/JMbayes/index.html',
                                                  class='lead'),
                                           tags$span(" to the world's largest active surveillance dataset called ", class='lead'),
                                           tags$a("PRIAS", target="_blank", href="http://www.prias-project.org/", class='lead'),
                                           tags$span(". The resulting prediction joint model is utilized in this web-application.", class='lead'),
                                           tags$hr(),
                                           tags$h2(tags$span("3. Cumulative-risk predictions are personalized and dynamic", class='label label-success')),
                                           tags$span("The risk predictions are ", class='lead'),
                                           tags$a("patient- and visit-specific [2]", href='https://doi.org/10.1002/bimj.201600238', target='_blank', class='lead'),
                                           tags$span(". That is, they are calculated separately for each patient, utilizing their personal history of PSA and biopsies.", class='lead'),
                                           tags$span("They are also updated on each follow-up visit with new data gathered until the follow-up.", class='lead'),
                                           tags$hr(),
                                           tags$h2(tags$span("4. Personalized biopsy schedules", class='label label-success')),
                                           tags$span("Personalized biopsy schedules utilize the cumulative-risk predictions of the patient over the entire follow-up period.", class='lead'),
                                           tags$span("More specifically, they schedule a biopsy at a follow-up visit if the cumulative-risk of Gleason Upgrading at a visit is higher than a certain threshold (e.g., 10% risk).", class='lead'),
                                           tags$span("Biopsies are scheduled sequentially, and after each decision of biopsy the cumulative-risk profile is updated to account for the possibility of not finding Gleason Upgrading during the biopsy.", class='lead'),
                                           tags$hr(),
                                           tags$h2(tags$span("5. Automatic Risk Threshold (Personalized) Schedule", class='label label-success')),
                                           tags$span("Along with fixed risk thresholds, this application also suggests an automatically chosen patient-specific risk threshold.", class='lead'), 
                                           tags$span("The automatic risk threshold minimizes the Euclidean distance between a risk-based schedule and an ideal schedule in a bi-dimensional space of the expected number of tests and expected time delay in detecting progression (figure below).", class='lead'),
                                           tags$span("An ideal schedule is a schedule that plans one biopsy exactly at the true time of Gleason upgrade.", class='lead'), 
                                           tags$span("The risk threshold with the least euclidean distance is chosen as the automatic risk threshold.", class='lead'), 
                                           tags$span("However, patients and doctors can vary this threshold and decide which resulting schedule suits their needs best by comparing the planned biopsies and expected time delay in detecting progression.", class='lead'),
                                           plotOutput("automatic_schedule_explanation", height = "600px"),
                                           tags$hr(),
                                           tags$h2(tags$span("6. Risk predictions and biopsies are planned only for a limited follow-up period", class='label label-success')),
                                           tags$span("We schedule biopsies only for a currently limited follow-up period in each cohort, because of limited study period of each cohort.", class='lead'), 
                                           tableOutput("table_max_pred_time"),
                                           tags$span("A compulsory biopsy was scheduled at the above maximum follow-up times in all schedules for meaningful comparison of their expected delays in detection of Gleason Upgrading.", class='lead'),
                                           tags$hr(),
                                           tags$h2(tags$span("7. Expected time delay in detecting Gleason Upgrading", class='label label-success')),
                                           tags$span("In active surveillance biopsies are conducted intermittently. Hence, Gleason Upgrading is always detected with a time delay as shown in the Figure below.", class='lead'),
                                           tags$span("For this patient biopsies were conducted at year 1 (Jan 2006), year 2 (Jan 2007) and year 4.5 (July 2009) since the start of active surveillance (Jan 2005).", class='lead'),
                                           tags$span("Since Gleason Upgrading occured at year 3.25 (Apr 2008), it could only be detected after a delay of 15 months at year 4.5 (Jul 2009).", class='lead'),
                                           fluidRow(column(10, plotOutput('delay_explanation_plot'))),
                                           tags$span("The smaller this delay, the larger is the window of opportunity for curative treatment.", class='lead'),
                                           tags$span("Using the cumulative-risk profile of a patient, for any given schedule of biopsies, we can estimate this time delay.", class='lead'),
                                           tags$span("In this web-application we calculate it for both personalized and fixed biopsy schedules.",class='lead'),
                                           tags$span("In this calculator, the expected time delay in detection of Gleason Upgrading is calculated for a scenario in which the patient obtains Gleason Upgrading before the date of the last planned biopsy in the schedule.", class='lead text-danger'),
                                           tags$hr(),
                                           tags$h2(tags$span("References", class='label label-primary')),
                                           tags$span("1. Rizopoulos, D. (2016). The R Package JMbayes for Fitting Joint Models for Longitudinal and Time-to-Event Data Using MCMC. Journal of Statistical Software, 72(7), 1 - 46. doi:http://dx.doi.org/10.18637/jss.v072.i07", class='lead'),
                                           tags$br(),
                                           tags$span("2. Rizopoulos, D., Molenberghs, G., & Lesaffre, E. M. (2017). Dynamic predictions with time-dependent covariates in survival analysis using joint modeling and landmarking. Biometrical Journal, 59(6), 1261-1276.", class='lead'))
                                  
                      )
                    )
                  )
)
)
