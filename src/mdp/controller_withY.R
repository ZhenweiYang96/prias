library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
#Source the method you want to use
source("src/mdp/tree_search/forward_search_with_Y.R")

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
  
  print(paste("Loaded file", file_num))
  
  fitted_JM = jointModelData$mvJoint_dre_psa_simDs
  testData = jointModelData$testData
  jointModelData$mvJoint_dre_psa_simDs = NULL
  
  LOWER_UPPER_PSA_LIMITS = by(data = testData$testDs$log2psaplus1,
                              INDICES = testData$testDs$visitTimeYears,
                              FUN = function(x){
                                q_one_two_third = quantile(x, probs = c(1/3,2/3))
                                ret = list(c(-Inf, q_one_two_third[1]),
                                           c(q_one_two_third[1], q_one_two_third[2]),
                                           c(q_one_two_third[2], Inf))
                              })
  names(LOWER_UPPER_PSA_LIMITS) = PSA_CHECK_UP_TIME[1:length(LOWER_UPPER_PSA_LIMITS)]
  
  sim_res = vector("list", length(time_to_biopsy_scales))
  names(sim_res) = names(time_to_biopsy_scales)
  for(i in 1:length(time_to_biopsy_scales)){
    print(paste("Time to biopsy scale", names(time_to_biopsy_scales)[i]))
    
    sim_res[[i]] = vector("list", length(max_depths))
    
    for(max_depth in max_depths){
      print(paste("Max depth:", max_depth))
      sim_res[[i]][[max_depth + 1]] = testData$testDs.id[,1:3]
      
      sim_res[[i]][[max_depth + 1]][,c('nb', 'offset')] = 
        foreach(pid=testData$testDs.id$P_ID, 
                .packages = c("splines", "JMbayes"), .combine = 'rbind') %dopar%{
                  source("src/mdp/common/simCommon.R")
                  debugSource("src/mdp/tree_search/forward_search_with_Y.R")
                  
                  patient_df = testData$testDs[testData$testDs$P_ID==pid,]
                  patient_df$psa_cat_data = F
                  progression_time = patient_df$progression_time[1]
                  time_to_biopsy_scale = time_to_biopsy_scales[i]
                  
                  set.seed(patient_df$P_ID[1])
                  
                  nb = 0
                  delay = -Inf
                  latest_biopsy_time = 0
                  decision_epoch = 1
                  while(decision_epoch<10 & delay<0){
                    pat_data = patient_df[patient_df$visitTimeYears <= decision_epoch,]
                    
                    DISCOUNT_FACTORS = DISCOUNT_FACTOR^((1:length(PSA_CHECK_UP_TIME))-nrow(pat_data))
                    names(DISCOUNT_FACTORS) = PSA_CHECK_UP_TIME
                    
                    act = selectAction(pat_data, current_decision_epoch = decision_epoch,
                                       latest_survival_time = latest_biopsy_time, earliest_failure_time = Inf,
                                       max_decision_epoch = min(10, decision_epoch + max_depth))
                    
                    if(act$optimal_action==BIOPSY){
                      latest_biopsy_time = decision_epoch
                      delay = decision_epoch - progression_time
                      decision_epoch = PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME >= (decision_epoch + 1)][1]
                      nb = nb + 1
                    }else{
                      decision_epoch = getNextDecisionEpoch(decision_epoch)
                    }
                  }
                  
                  return(c('nb'=nb, 'offset'=delay))
                }
      
      save(sim_res, file = paste0("Rdata/mdp/schedule_by_time_Y/sim_res_",file_num,"_", 
                                  names(time_to_biopsy_scales)[i], ".Rdata"))
    }
  }
  rm(fitted_JM)
}

stopCluster(ct)
