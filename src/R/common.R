library(ggplot2)
library(coda)
library(ggmcmc)
library(parallel)
library(doParallel)
library(JMbayes)

ticksX = function(from=0, max, by, labels=waiver()){
  scale_x_continuous(breaks = seq(from, max, by = by), labels = labels)
}

ticksY = function(from=0, max, by, labels = waiver()){
  scale_y_continuous(breaks = seq(from, max, by = by), labels=waiver())
}

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#Chapter 5, section 5.2
#shiny::runGitHub("Repeated_Measurements", "drizopoulos")

