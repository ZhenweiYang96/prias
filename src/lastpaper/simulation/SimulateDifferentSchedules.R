rm(list=ls())

#One choice of a seed provides one round of simulation
seed = 2021

#Number of parallel cores to use
max_cores = 4

source("PrepareSimStudyJointModel.R")

#Choose the ID of the patient for whom we intend to conduct schedules
#testId range between 751 and 1000 (both included)
testId=751

testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID == testId,]
testDs.id = jointModelData$testData$testDs.id[jointModelData$testData$testDs.id$P_ID == testId,]
progression_time = testDs.id$progression_time

ANNUAL = "Annual"
PRIAS = "PRIAS"

KAPPApt05 = "Risk: 5%"
KAPPApt1 = "Risk: 10%"
KAPPApt15 = "Risk: 15%"

#M is the number of MCMC simulations
M=500
FIRST_DECISION_VISIT_NR = 5
MIN_BIOPSY_GAP = 1

set.seed(seed)

print(paste("**** Patient ID: ", testId, 
            " with progression time: ", progression_time,
            " from simulation with seed: ", seed))

#special case
if(progression_time<=1){
  print("progression time is less than 1 year")
  annual_biopsies = prias_biopsies = 1
  print('done annual schedule and prias schedule')
}else{
  #annual and biennial biopsies
  annual_biopsies = seq(1, ceiling(progression_time), by = 1)
  
  print('done annual schedule')
  
  #PRIAS schedule
  prias_biopsies = c(1)
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifPRIASBiopsy(patient_df, cur_visit_time, max(prias_biopsies,0))){
      prias_biopsies = c(prias_biopsies, cur_visit_time)
    }
    
    if(max(prias_biopsies, 0) >= progression_time){
      break
    }
  }
  #check if progression could not be detected within 10 years despite being progressed
  if(max(prias_biopsies, 0) < progression_time){
    prias_biopsies = c(prias_biopsies, MAX_FAIL_TIME)
  }
  print('done PRIAS schedule')
}

#Risk based biopsies
fixed_risk_biopsies = lapply(c(0.05, 0.1, 0.15), function(threshold){
  source("prediction_psa_dre.R")
  risk_biopsies = c()
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifFixedRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                              patient_data = patient_df, 
                              cur_visit_time = cur_visit_time, 
                              last_biopsy_time = max(risk_biopsies, 0),
                              min_biopsy_gap = MIN_BIOPSY_GAP,
                              threshold = threshold,
                              M = M)){
      risk_biopsies = c(risk_biopsies, cur_visit_time)
    }
    
    if(max(risk_biopsies, 0) >= progression_time){
      break
    }
  }
  if(max(risk_biopsies, 0) < progression_time){
    risk_biopsies = c(risk_biopsies, MAX_FAIL_TIME)
  }
  
  print(paste0('done fixed risk biopsy with threshold: ', threshold*100, "%"))
  return(risk_biopsies)
})

#automatically selected risk threshold
auto_risk_biopsies = lapply(c(0.75, Inf), function(delay_limit){
  source("prediction_psa_dre_randEff_reuse.R")
  automatic_kappa_biopsies = c()
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifAutomaticRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                                  patient_data = patient_df, 
                                  cur_visit_time = cur_visit_time, 
                                  last_biopsy_time = max(automatic_kappa_biopsies, 0),
                                  min_biopsy_gap = MIN_BIOPSY_GAP, delay_limit = delay_limit,
                                  M = M, horizon = MAX_FAIL_TIME, use_exact_delay = T)){
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

KAPPA_auto_pt75 = "Risk: Auto (0.75)"
KAPPA_auto_Inf = "Risk: Auto (Inf)"

scheduleNames = c(ANNUAL, BIENNIAL, PRIAS, 
                  KAPPApt05, KAPPApt1, KAPPApt15, 
                  KAPPA_auto_pt75, KAPPA_auto_Inf)

scheduleLengths = c(length(annual_biopsies), length(biennial_biopsies), length(prias_biopsies),
                    sapply(fixed_risk_biopsies, length), sapply(auto_risk_biopsies, length))

biopsy_time = c(annual_biopsies, biennial_biopsies, prias_biopsies,
                do.call('c', fixed_risk_biopsies),
                do.call('c', auto_risk_biopsies))

#Now we combine all of these into a data frame for this patient
biopsyDf = data.frame(seed = seed, P_ID = testId, 
                      progression_time = progression_time,
                      progressed = testDs.id$progressed,
                      schedule = rep(scheduleNames,scheduleLengths),
                      biopsy_time = biopsy_time)