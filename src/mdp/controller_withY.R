library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction_psa_cat.R")
source("src/mdp/tree_search/forward_search_with_Y.R")

#Decide the number of cores
max_cores = 3
N_MCMC_ITER = 0

#For loop parameters
nDs = 10
max_depths = c(0,1, 4,3,2)
discount_factors = 1
time_to_biopsy_scales = unique(c(10/0.5, 7/0.7, 5/1.1, 4/1.8, 4/0.7, 3/0.8))
names(time_to_biopsy_scales) = c("Annual", "1pt5years", "Biennial", "Triennial",
                                 "PRIAS_risk10pc", "risk15pc")

ct = makeCluster(max_cores, type = "FORK")
registerDoParallel(ct)
for(file_num in nDs:1){
  file = list.files(path='/home/a_tomer/Data/final_res_2nd_paper/', full.names = T)[file_num]
  #Load a simulation
  load(file)
  
  print(paste("Loaded file:", file_num))
  
  fitted_JM = jointModelData$mvJoint_dre_psa_simDs
  testData = jointModelData$testData
  jointModelData$mvJoint_dre_psa_simDs = NULL
  
  LOWER_UPPER_PSA_LIMITS = by(data = testData$testDs$log2psaplus1,
                              INDICES = testData$testDs$visitTimeYears,
                              FUN = function(x){
                                qtiles = quantile(x, probs = seq(0, 1, 1/NR_DISCRETIZED_PSA))
                                qtiles[1] = -Inf
                                qtiles[length(qtiles)] = Inf
                                lapply(1:NR_DISCRETIZED_PSA, function(x){
                                  c(qtiles[x], qtiles[x+1]) 
                                })
                              })
  names(LOWER_UPPER_PSA_LIMITS) = PSA_CHECK_UP_TIME[1:length(LOWER_UPPER_PSA_LIMITS)]
  
  for(DISCOUNT_FACTOR in discount_factors){
    print(paste("Discount factor:", DISCOUNT_FACTOR))
    
    sim_res = vector("list", length(time_to_biopsy_scales))
    names(sim_res) = names(time_to_biopsy_scales)
    for(i in 1:length(time_to_biopsy_scales)){
      print(paste("Time to biopsy scale", names(time_to_biopsy_scales)[i]))
      
      sim_res[[i]] = vector("list", length(max_depths))
      
      for(max_depth in max_depths){
        print(paste("Max depth:", max_depth))
        sim_res[[i]][[max_depth + 1]] = testData$testDs.id[, c("P_ID", "Age", "progression_time")]
        
        sim_res[[i]][[max_depth + 1]][, c('nb', 'offset')] = 
          foreach(pid=testData$testDs.id$P_ID, 
                  .packages = c("splines", "JMbayes"), .combine = 'rbind') %do%{
                    
                    patient_df = testData$testDs[testData$testDs$P_ID==pid,]
                    patient_df$psa_cat_data = F
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
                      
                      max_decision_epoch = min(MAX_FOLLOW_UP_TIME, decision_epoch + max_depth)
                      act = selectAction(pat_data, current_decision_epoch = decision_epoch,
                                         latest_survival_time = latest_biopsy_time, earliest_failure_time = Inf,
                                         max_decision_epoch = max_decision_epoch)
                      
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
        
        save(sim_res, file = paste0("Rdata/mdp/schedule_by_time_Y/sim_res_",file_num,"_", 
                                    names(time_to_biopsy_scales)[i], ".Rdata"))
      }
    }
  }
  rm(fitted_JM)
}

stopCluster(ct)
