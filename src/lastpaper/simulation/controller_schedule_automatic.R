args = commandArgs(trailingOnly=TRUE)

seed = as.numeric(args[1])
testId = as.numeric(args[2])

library(JMbayes)
library(splines)
library(survival)
library(MASS)

load(paste0("Rdata/lastpaper/simulation/light/jointModelData_seed_", seed, "_t3.Rdata"))
source("src/lastpaper/prediction_psa_dre_randEff_reuse.R")
source("src/lastpaper/simulation/scheduleCreator.R")

KAPPA_auto_pt5 = "Risk: Auto (0.5)"
KAPPA_auto_1 = "Risk: Auto (1)"
KAPPA_auto_Inf = "Risk: Auto (Inf)"

M=500
MAX_FAIL_TIME = 10
FIRST_DECISION_VISIT_NR = 5
MIN_BIOPSY_GAP = 1

testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID == testId,]
testDs.id = jointModelData$testData$testDs.id[jointModelData$testData$testDs.id$P_ID == testId,]
progression_time = testDs.id$progression_time

set.seed(seed)

print(paste("**** Patient ID: ", testId, 
            " with progression time: ", progression_time,
            " from simulation with seed: ", seed))

#automatically selected risk threshold
auto_risk_biopsies = lapply(c(0.5, 1, Inf), function(delay_limit){
  automatic_kappa_biopsies = c()
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifAutomaticRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                                  patient_data = patient_df, 
                                  cur_visit_time = cur_visit_time, 
                                  last_biopsy_time = max(automatic_kappa_biopsies, 0),
                                  min_biopsy_gap = MIN_BIOPSY_GAP, delay_limit = delay_limit,
                                  M = M, horizon = MAX_FAIL_TIME, use_restricted_delay = T)){
      automatic_kappa_biopsies = c(automatic_kappa_biopsies, cur_visit_time)
    }
    
    if(max(automatic_kappa_biopsies, 0) >= progression_time){
      break
    }
  }
  if(max(automatic_kappa_biopsies, 0) < progression_time){
    automatic_kappa_biopsies = c(automatic_kappa_biopsies, MAX_FAIL_TIME)
  }
  print(paste0('done automatic risk biopsy with delay: ', delay_limit, " years"))
  return(automatic_kappa_biopsies)
})

scheduleNames = c(KAPPA_auto_pt5, KAPPA_auto_1, KAPPA_auto_Inf)
scheduleLengths = sapply(auto_risk_biopsies, length)

biopsy_time = do.call('c', auto_risk_biopsies)
#Now we combine all of these into a data frame for this patient
biopsyDf = data.frame(seed = seed, P_ID = testId, 
                      progression_time = progression_time,
                      progressed = testDs.id$progressed,
                      schedule = rep(scheduleNames,scheduleLengths),
                      biopsy_time = biopsy_time)

print("Saving now")
patient_file_name = paste0("auto_seed_", seed, "_P_ID_", testId)
save(biopsyDf, file=paste0("Rdata/lastpaper/simulation/results/", patient_file_name, ".Rdata"))
