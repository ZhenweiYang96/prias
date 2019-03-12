library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction.R")
source("src/mdp/tree_search/forward_search_no_Y.R")

#For DESPOT set these two
#N_DESPOT_SCENARIOS = 100
#DESPOT_TREE = list()

max_cores = 3

max_depths = 0:4
time_to_biopsy_scales = unique(c(10/0.5, 7/0.7, 5/1.1, 4/1.8, 4/0.7, 3/0.8))
names(time_to_biopsy_scales) = c("Annual", "1pt5years", "Biennial", "Triennial",
                                 "PRIAS_risk10pc", "risk15pc")

ct = makeCluster(max_cores, type = "FORK")
registerDoParallel(ct)

nDs = 10
for(file_num in nDs:1){
  
  file = list.files(path='/home/a_tomer/Data/final_res_2nd_paper/', full.names = T)[file_num]
  #Load a simulation
  load(file)
  
  fitted_JM = jointModelData$mvJoint_dre_psa_simDs
  jointModelData$mvJoint_dre_psa_simDs = NULL
  
  sim_res = vector("list", length(time_to_biopsy_scales))
  names(sim_res) = names(time_to_biopsy_scales)
  for(i in 1:length(time_to_biopsy_scales)){
    sim_res[[i]] = vector("list", length(max_depths))
    
    for(max_depth in max_depths){
      sim_res[[i]][[max_depth + 1]] = jointModelData$testData$testDs.id[,1:3]
      
      sim_res[[i]][[max_depth + 1]][,c('nb', 'offset')] = 
        foreach(pid=jointModelData$testData$testDs.id$P_ID, 
                .packages = c("splines", "JMbayes"), .combine = 'rbind') %dopar%{
                  
                  patient_df = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID==pid,]
                  progression_time = patient_df$progression_time[1]
                  time_to_biopsy_scale = time_to_biopsy_scales[i]
                  
                  set.seed(patient_df$P_ID[1])
                  
                  nb = 0
                  delay = -Inf
                  latest_biopsy_time = 0
                  decision_epoch = 1
                  while(decision_epoch<MAX_FOLLOW_UP_TIME & delay<0){
                    pat_data = patient_df[patient_df$visitTimeYears <= decision_epoch,]
                    
                    DISCOUNT_FACTORS = DISCOUNT_FACTOR^((1:length(BIOPSY_TEST_TIMES))-nrow(pat_data))
                    names(DISCOUNT_FACTORS) = BIOPSY_TEST_TIMES
                    
                    act = selectAction(pat_data, current_decision_epoch = decision_epoch,
                                       latest_survival_time = latest_biopsy_time, earliest_failure_time = Inf,
                                       max_decision_epoch = min(MAX_FOLLOW_UP_TIME, decision_epoch + max_depth))
                    
                    if(act$optimal_action==BIOPSY){
                      latest_biopsy_time = decision_epoch
                      delay = decision_epoch - progression_time
                      decision_epoch = BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES >= (decision_epoch + MIN_BIOPSY_GAP)][1]
                      nb = nb + 1
                    }else{
                      decision_epoch = getNextDecisionEpoch(decision_epoch)
                    }
                  }
                  
                  return(c('nb'=nb, 'offset'=delay))
                }
      
      save(sim_res, file = paste0("Rdata/mdp/schedule_by_time/sim_res_",file_num,"_", 
                                  names(time_to_biopsy_scales)[i], ".Rdata"))
    }
  }
  rm(fitted_JM)
}

stopCluster(ct)
