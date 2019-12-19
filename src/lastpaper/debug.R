library(JMbayes)
library(splines)
library(survival)
library(MASS)

#View(biopsyDf_summary[biopsyDf_summary$progression_time==10 & 
#                        biopsyDf_summary$schedule=="Risk: Auto (1)",])

seed = 2001
testId = 890

load(paste0("Rdata/lastpaper/simulation/light/jointModelData_seed_", seed, "_t3.Rdata"))
source("src/lastpaper/prediction_psa_dre.R")
source("src/lastpaper/scheduleCreator.R")

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
#progression_time = 10
delay_limit = 0.5
tree_biopsies = c()
for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
  patient_df = testDs[1:row_num,]
  cur_visit_time = testDs$year_visit[row_num]
  
  if(ifAutomaticFullTreeBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                                patient_data = patient_df, 
                                cur_visit_time = cur_visit_time, 
                                last_biopsy_time = max(tree_biopsies, 0),
                                min_biopsy_gap = MIN_BIOPSY_GAP, delay_limit = delay_limit,
                                M = M, horizon = MAX_FAIL_TIME)){
    tree_biopsies = c(tree_biopsies, cur_visit_time)
  }
  
  if(max(tree_biopsies,0) >= progression_time){
    break
  }
}
print(tree_biopsies)