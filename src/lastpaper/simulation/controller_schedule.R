args = commandArgs(trailingOnly=TRUE)

seed = as.numeric(args[1])
testId = as.numeric(args[2])

library(JMbayes)
library(splines)
library(survival)
library(MASS)

load(paste0("Rdata/lastpaper/simulation/light/jointModelData_seed_", seed, "_t3.Rdata"))
source("src/lastpaper/prediction_psa_dre.R")
source("src/lastpaper/scheduleCreator.R")

ANNUAL = "Annual"
BIENNIAL = "Biennial"
PRIAS = "PRIAS"

KAPPApt05 = "Risk: 5%"
KAPPApt1 = "Risk: 10%"
KAPPApt15 = "Risk: 15%"
KAPPAautomatic_pt5 = "Risk: Auto:"

methodNames = c(ANNUAL, BIENNIAL, PRIAS, KAPPApt05, KAPPApt1, KAPPApt15, KAPPAautomatic)

M=500
MAX_FAIL_TIME = 10
FIRST_DECISION_VISIT_NR = 9
MIN_BIOPSY_GAP = 1

print(paste("**** Doing for patient ID ", testId, 
            " from simulation with seed ", seed))

testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID == testId,]
testDs.id = jointModelData$testData$testDs.id[jointModelData$testData$testDs.id$P_ID == testId,]

#special case
if(testDs.id$progression_time<=1){
  
}

#annual and biennial biopsies
annual_biopsies = seq(1, ceiling(testDs.id$progression_time), by = 1)
biennial_biopsies = c(1,3,5,7,9,10)
biennial_biopsies = biennial_biopsies[1:(which(biennial_biopsies-testDs.id$progression_time >= 0)[1])]

#PRIAS schedule
prias_biopsies = c(1)
for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
  patient_df = testDs[1:row_num,]
  cur_visit_time = testDs$year_visit[row_num]
  
  if(ifPRIASBiopsy(patient_df, cur_visit_time, tail(prias_biopsies,1))){
    prias_biopsies = c(prias_biopsies, cur_visit_time)
  }
  
  if(tail(prias_biopsies,1) >= testDs.id$progression_time){
    break
  }
}

#Risk based biopsies
fixed_risk_biosies = lapply(c(0.05, 0.1, 0.15), function(threshold){
  risk_biopsies = c(1)
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifFixedRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                              patient_data = patient_df, 
                              cur_visit_time = cur_visit_time, 
                              last_biopsy_time = tail(risk_biopsies,1),
                              min_biopsy_gap = MIN_BIOPSY_GAP,
                              threshold = threshold,
                              M = M)){
      risk_biopsies = c(risk_biopsies, cur_visit_time)
    }
    
    if(tail(risk_biopsies,1) >= testDs.id$progression_time){
      break
    }
  }
  return(risk_biopsies)
})
names(fixed_risk_biosies) = c(KAPPApt05, KAPPApt1, KAPPApt15)

#Use different delays...better than being sorry later
c(0.5, 1, 2, Inf)

#automatically selected risk threshold
automatic_kappa_biopsies = c(1)
for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
  patient_df = testDs[1:row_num,]
  cur_visit_time = testDs$year_visit[row_num]
  
  if(ifAutomaticRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                                patient_data = patient_df, 
                                cur_visit_time = cur_visit_time, 
                                last_biopsy_time = tail(automatic_kappa_biopsies,1),
                                min_biopsy_gap = MIN_BIOPSY_GAP, delay_limit = 0.5,
                                M = M, horizon = MAX_FAIL_TIME)){
    automatic_kappa_biopsies = c(automatic_kappa_biopsies, cur_visit_time)
  }
  
  if(tail(automatic_kappa_biopsies,1) >= testDs.id$progression_time){
    break
  }
}
