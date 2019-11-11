args = commandArgs(trailingOnly = F)

sim_seed = as.numeric(args[1])

#9 corresponds to year 2, and 5 corresponds to year 1
FIRST_DECISION_VISIT_NR = 9
MIN_BIOPSY_GAP = 1

M = 750

library(JMbayes)
library(splines)

print(paste0("Loading simulation with seed ", sim_seed))
load(paste0("Rdata/lastpaper/sims/sim_seed_", sim_seed, ".Rdata"))

source('src/lastpaper/minDistThreshold.R')
source('src/lastpaper/prediction_only_psa.R')
source('src/lastpaper/scheduleCreatorCacheBased.R')

nb1_to_delay = as.numeric(args[2])
delay_limit = as.numeric(args[3])
MAX_FOLLOW_UP = as.numeric(args[4])

print(paste0("Maximum horizon is ", MAX_FOLLOW_UP))
print(paste0("Maximum limit on delay is ", delay_limit))
print(paste0("1 biopsy is equal to ", nb1_to_delay, " units of delay"))

testDs = sim_res$testData$testDs
testDs = testDs[testDs$year_visit <= MAX_FOLLOW_UP,]
testDs.id = sim_res$testData$testDs.id

timesPerSubject = max(testDs$visitNumber)

nbCol = paste0("nb_", nb1_to_delay)
delayCol = paste0("delay_", nb1_to_delay)

testDs.id[, nbCol] = 0
testDs.id[, delayCol] = NA

set.seed(sim_seed)
for(testId in testDs.id$P_ID){
  
  print(paste("Doing for patient:", testId))
  
  curVisitNr = FIRST_DECISION_VISIT_NR
  
  patient_data = testDs[testDs$P_ID == testId,]
  patient_data$gleason_sum[1] = 6
  
  idFilter = testDs.id$P_ID == testId
  progression_time = testDs.id$progression_time[idFilter]
  
  latest_survival_time = 0
  
  while(curVisitNr <= timesPerSubject){
    cur_visit_time = patient_data$year_visit[curVisitNr]
    if(cur_visit_time - latest_survival_time >= MIN_BIOPSY_GAP){
      
      decision = minDistScheduleDecision(object = sim_res$mvJoint_psa_simDs, 
                                         patient_data = patient_data[1:curVisitNr,], 
                                         cur_visit_time = cur_visit_time,
                                         latest_survival_time = latest_survival_time,
                                         delay_limit = delay_limit,
                                         nb1_to_delay = nb1_to_delay,
                                         weight_by_horizon_risk = F)
      
      if(decision$result == T){
        patient_data$gleason_sum[curVisitNr] = 6
        testDs.id[idFilter, nbCol] = testDs.id[idFilter, nbCol] + 1
        testDs.id[idFilter, delayCol] = cur_visit_time - progression_time
        
        latest_survival_time = cur_visit_time
        
        if(testDs.id[idFilter, delayCol] >= 0){
          break
        }
      }
    }
    
    curVisitNr = curVisitNr + 1
  }
}

sim_res$testData$testDs.id = testDs.id

save(sim_res, file=paste("Rdata/lastpaper/schedule_res", sim_seed, nb1_to_delay, delay_limit, MAX_FOLLOW_UP, ".Rdata", sep = "_"))