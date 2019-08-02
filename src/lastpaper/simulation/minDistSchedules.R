rm(list = ls())
library(JMbayes)
library(splines)
load("Rdata/lastpaper/sims/sim_seed_2019.Rdata")
source("src/lastpaper/prediction.R")
source("src/lastpaper/minDistThreshold.R")
source("src/clinical_gap3/scheduleCreator.R")

timesPerSubject = max(sim_res$testData$testDs$visitNumber)

set.seed(2019)
for(nb1_to_delay in seq(1.75, 2, 0.25)){
  nbCol = paste0("nb_", nb1_to_delay)
  delayCol = paste0("delay_", nb1_to_delay)
  
  sim_res$testData$testDs.id[, nbCol] = 0
  sim_res$testData$testDs.id[, delayCol] = NA
  
  for(testId in sim_res$testData$testDs.id$P_ID){
    
    print(paste("Doing for patient:", testId))
    curVisitNr = 5
    patient_data = sim_res$testData$testDs[sim_res$testData$testDs$P_ID == testId,]
    patient_data$gleason_sum[1] = 6
    
    idFilter = sim_res$testData$testDs.id$P_ID == testId
    progression_time = sim_res$testData$testDs.id$progression_time[idFilter]
    
    latest_survival_time = 0
    
    while(curVisitNr <= timesPerSubject){
      cur_visit_time = patient_data$year_visit[curVisitNr]
      if(cur_visit_time - latest_survival_time >= 1){
        decision = minDistScheduleDecision(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], 
                                           cur_visit_time = cur_visit_time,
                                           latest_survival_time = latest_survival_time,
                                           nb1_to_delay=nb1_to_delay,
                                           weight_by_horizon_risk = F)
        
        if(decision==T){
          patient_data$gleason_sum[curVisitNr] = 6
          sim_res$testData$testDs.id[idFilter, nbCol] = sim_res$testData$testDs.id[idFilter, nbCol] + 1
          sim_res$testData$testDs.id[idFilter, delayCol] = cur_visit_time - progression_time
          latest_survival_time = cur_visit_time
          if(sim_res$testData$testDs.id[idFilter, delayCol]>=0){
            break;
          }
        }
      }
      
      curVisitNr = curVisitNr + 1
    }
  }
  
  save.image("Rdata/lastpaper/sim_seed_2019_sched_wtbyhorizonriskFALSE.Rdata")
}
