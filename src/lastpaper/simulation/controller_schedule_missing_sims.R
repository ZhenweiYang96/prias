args = commandArgs(trailingOnly=TRUE)

load("Rdata/lastpaper/simulation/missing_sims.Rdata")

patient = missing_sims[as.numeric(args[1]),]

seed = patient$seed
testId = patient$P_ID
schedule = as.character(patient$schedule)

library(JMbayes)
library(splines)
library(survival)
library(MASS)

fixed_threshold = c(0.05,0.1,0.15)
names(fixed_threshold) = c("Risk: 5%", "Risk: 10%", "Risk: 15%")

automatic_delay_limit = c(0.75, Inf)
names(automatic_delay_limit) = c("Risk: Auto (0.75)", "Risk: Auto (Inf)")


load(paste0("Rdata/lastpaper/simulation/light/jointModelData_seed_", seed, "_t3.Rdata"))

if(schedule %in% names(automatic_delay_limit)){
  source("src/lastpaper/prediction_psa_dre_randEff_reuse.R")
}else{
  source("src/lastpaper/prediction_psa_dre.R")  
}
source("src/lastpaper/simulation/scheduleCreator.R")

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
            " from simulation with seed: ", seed,
            " schedule: ", schedule))

#special case
if(progression_time<=1){
  annual_biopsies = biennial_biopsies = prias_biopsies = 1
}else{
  #annual and biennial biopsies
  annual_biopsies = seq(1, ceiling(progression_time), by = 1)
  biennial_biopsies = c(1,3,5,7,9,10)
  biennial_biopsies = biennial_biopsies[1:(which(biennial_biopsies-progression_time >= 0)[1])]
  
  if(schedule=="PRIAS"){
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
  }
}

risk_biopsies = c()
if(schedule %in% names(fixed_threshold)){
  threshold = fixed_threshold[schedule]
  
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
}

automatic_kappa_biopsies = c()
if(schedule %in% names(automatic_delay_limit)){
  delay_limit = automatic_delay_limit[schedule]
  
  for(row_num in FIRST_DECISION_VISIT_NR:nrow(testDs)){
    patient_df = testDs[1:row_num,]
    cur_visit_time = testDs$year_visit[row_num]
    
    if(ifAutomaticRiskBasedBiopsy(object = jointModelData$mvJoint_dre_psa_simDs,
                                  patient_data = patient_df, 
                                  cur_visit_time = cur_visit_time, 
                                  last_biopsy_time = max(automatic_kappa_biopsies, 0),
                                  min_biopsy_gap = MIN_BIOPSY_GAP, delay_limit = delay_limit,
                                  M = M, horizon = MAX_FAIL_TIME)){
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
}

switch(schedule, 
       "Annual"={
         biopsy_time = annual_biopsies
       },
       "Biennial"={
         biopsy_time = biennial_biopsies
       },
       "PRIAS"={
         biopsy_time = prias_biopsies
       },
       "Risk: 5%"={
         biopsy_time = risk_biopsies
       },
       "Risk: 10%"={
         biopsy_time = risk_biopsies
       },
       "Risk: 15%"={
         biopsy_time = risk_biopsies
       },
       "Risk: Auto (0.75)"={
         biopsy_time = automatic_kappa_biopsies
       },
       "Risk: Auto (Inf)"={
         biopsy_time = automatic_kappa_biopsies
       }
)

biopsyDf = data.frame(seed = seed, P_ID = testId, 
                      progression_time = progression_time,
                      progressed = testDs.id$progressed,
                      schedule = schedule,
                      biopsy_time = biopsy_time)

print("Saving now")
patient_file_name = paste0("seed_", seed, "_P_ID_", testId, "_schedule_", schedule)
save(biopsyDf, file=paste0("Rdata/lastpaper/simulation/results/", patient_file_name, ".Rdata"))
