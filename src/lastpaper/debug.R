seed = 2021
testId = 751

library(JMbayes)
library(splines)
library(survival)
library(MASS)

load(paste0("Rdata/lastpaper/simulation/light/jointModelData_seed_", seed, "_t3.Rdata"))
source("src/lastpaper/prediction_psa_dre.R")
#source("src/lastpaper/scheduleCreator.R")

M=500
MAX_FAIL_TIME = 10
FIRST_DECISION_VISIT_NR = 9
MIN_BIOPSY_GAP = 1

testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID == testId,]
testDs.id = jointModelData$testData$testDs.id[jointModelData$testData$testDs.id$P_ID == testId,]
progression_time = testDs.id$progression_time

set.seed(seed)

print(paste("**** Patient ID: ", testId, 
            " with progression time: ", progression_time,
            " from simulation with seed: ", seed))


#automatically selected risk threshold
delay = Inf
automatic_kappa_biopsies = c(1)
for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
  patient_df = testDs[1:row_num,]
  cur_visit_time = testDs$year_visit[row_num]
  
  if(ifAutomaticRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                                patient_data = patient_df, 
                                cur_visit_time = cur_visit_time, 
                                last_biopsy_time = tail(automatic_kappa_biopsies,1),
                                min_biopsy_gap = MIN_BIOPSY_GAP, delay_limit = delay,
                                M = M, horizon = MAX_FAIL_TIME)){
    automatic_kappa_biopsies = c(automatic_kappa_biopsies, cur_visit_time)
  }
  
  if(tail(automatic_kappa_biopsies,1) >= progression_time){
    break
  }
}
if(tail(automatic_kappa_biopsies,1) < progression_time){
  automatic_kappa_biopsies = c(automatic_kappa_biopsies, MAX_FAIL_TIME)
}
print(paste0('done automatic risk biopsy with delay: ', delay, " years"))
return(automatic_kappa_biopsies)

