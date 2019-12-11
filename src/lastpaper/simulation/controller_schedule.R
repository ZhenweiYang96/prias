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

KAPPA_auto_pt5 = "Risk: Auto (0.5)"
KAPPA_auto_1 = "Risk: Auto (1)"
KAPPA_auto_2 = "Risk: Auto (2)"
KAPPA_auto_Inf = "Risk: Auto (Inf)"

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

#special case
if(progression_time<=1){
  print("progression time is less than 1 year")
  annual_biopsies = biennial_biopsies = prias_biopsies = 1
  fixed_risk_biopsies = list(1,1,1)
  auto_risk_biopsies = list(1,1,1,1)
}else{
  #annual and biennial biopsies
  annual_biopsies = seq(1, ceiling(progression_time), by = 1)
  biennial_biopsies = c(1,3,5,7,9,10)
  biennial_biopsies = biennial_biopsies[1:(which(biennial_biopsies-progression_time >= 0)[1])]
  
  print('done annual and biennial')
  
  #PRIAS schedule
  prias_biopsies = c(1)
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifPRIASBiopsy(patient_df, cur_visit_time, tail(prias_biopsies,1))){
      prias_biopsies = c(prias_biopsies, cur_visit_time)
    }
    
    if(tail(prias_biopsies,1) >= progression_time){
      break
    }
  }
  #check if progression could not be detected within 10 years despite being progressed
  if(tail(prias_biopsies,1) < progression_time){
    prias_biopsies = c(prias_biopsies, MAX_FAIL_TIME)
  }
  print('done PRIAS')
  
  #Risk based biopsies
  fixed_risk_biopsies = lapply(c(0.05, 0.1, 0.15), function(threshold){
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
      
      if(tail(risk_biopsies,1) >= progression_time){
        break
      }
    }
    if(tail(risk_biopsies,1) < progression_time){
      risk_biopsies = c(risk_biopsies, MAX_FAIL_TIME)
    }
    
    print(paste0('done fixed risk biopsy with threshold: ', threshold*100, "%"))
    return(risk_biopsies)
  })
  
  #automatically selected risk threshold
  auto_risk_biopsies = lapply(c(0.5, 1, 2, Inf), function(delay){
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
  })
}

scheduleNames = c(ANNUAL, BIENNIAL, PRIAS, 
                  KAPPApt05, KAPPApt1, KAPPApt15, 
                  KAPPA_auto_pt5, KAPPA_auto_1, KAPPA_auto_2, KAPPA_auto_Inf)
scheduleLengths = c(length(annual_biopsies), length(biennial_biopsies), length(prias_biopsies),
                    sapply(fixed_risk_biopsies, length), sapply(auto_risk_biopsies, length))
biopsy_time = c(annual_biopsies, biennial_biopsies, prias_biopsies,
                do.call('c', fixed_risk_biopsies), do.call('c', auto_risk_biopsies))

#Now we combine all of these into a data frame for this patient
biopsyDf = data.frame(seed = seed, P_ID = testId, 
                      progression_time = progression_time,
                      progressed = testDs.id$progressed,
                      schedule = rep(scheduleNames,scheduleLengths),
                      biopsy_time = biopsy_time)

print("Saving now")
patient_file_name = paste0("seed_", seed, "_P_ID_", testId)
save(biopsyDf, file=paste0("Rdata/lastpaper/simulation/results/", patient_file_name, ".Rdata"))
