source("src/lastpaper/simulation/simCommon_psa.R")

nSub = 1200
bNames = paste0("b", 0:4)
set.seed(2019)
sim_data= generateSimulationData(nSub, bNames = bNames)
sim_res = fitJointModelOnNewData(sim_data$simDs, sim_data$simDs.id, nSub * 0.75)

timesPerSubject = max(sim_res$testData$testDs$visitNumber)
sim_res$testData$testDs$gleason_sum = NA
sim_res$testData$testDs$psa = 2^sim_res$testData$testDs$log2psaplus1 - 1

set.seed(2019)
for(nb_offset_ratio in seq(0, 2, 0.25)){
  sim_res$testData$testDs.id$nb_mindist = 0
  sim_res$testData$testDs.id$delay_mindist = NA
  
  nbCol = paste0("nb_nor", nb_offset_ratio)
  delayCol = paste0("delay_nor", nb_offset_ratio)
  
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
                                           nb_offset_ratio=nb_offset_ratio,
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
  
  save.image("Rdata/lastpaper/sim_seed_2019_sched_weight.Rdata")
}

#Risk based schedule
for(threshold in c(0.05, 0.1, 0.15)){
  
  nbCol = paste0("nb_", threshold)
  delayCol = paste0("delay_", threshold)
  
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
        curCumRisk = mean(1-getExpectedFutureOutcomes(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], latest_survival_time, Inf,
                                                      survival_predict_times = cur_visit_time)$predicted_surv_prob)
        decision = curCumRisk >= threshold
        
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
  
  save.image("Rdata/lastpaper/sim_seed_2019_sched.Rdata")
}

