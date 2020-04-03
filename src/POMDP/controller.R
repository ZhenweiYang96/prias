library(JMbayes)
library(survival)
library(splines)
library(ggplot2)
library(doParallel)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction.R")
source("src/mdp/tree_search/forward_search_no_Y.R")

max_cores = 3
N_MCMC_ITER = 300
N_DESPOT_SCENARIOS = 500
discount_factors = c(0.5, 1)
max_biopsies_vec = c(Inf)
max_depths = c(0, 1, 3, 4)

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)
reward_matrix = as.matrix(expand.grid(1,
                                      -c(3/12, 6/12), 
                                      0,
                                      -c(0/12, 6/12)))
colnames(reward_matrix) = reward_names

ct = makeCluster(max_cores)
registerDoParallel(ct)

t1 = Sys.time()
nDs = 1
for(file_num in nDs:1){
  
  file = list.files(path='/home/a_tomer/Data/final_res_2nd_paper/', full.names = T)[file_num]
  load(file)
  
  print(paste("Loaded file:", file_num))
  
  fitted_JM = jointModelData$mvJoint_dre_psa_simDs
  testData = jointModelData$testData
  jointModelData$mvJoint_dre_psa_simDs = NULL
  
  # LOWER_UPPER_PSA_LIMITS = by(data = testData$testDs$log2psaplus1,
  #                             INDICES = testData$testDs$visitTimeYears,
  #                             FUN = function(x){
  #                               qtiles = quantile(x, probs = seq(0, 1, 1/NR_DISCRETIZED_PSA))
  #                               qtiles[1] = -Inf
  #                               qtiles[length(qtiles)] = Inf
  #                               lapply(1:NR_DISCRETIZED_PSA, function(x){
  #                                 c(qtiles[x], qtiles[x+1]) 
  #                               })
  #                             })
  # names(LOWER_UPPER_PSA_LIMITS) = PSA_CHECK_UP_TIME[1:length(LOWER_UPPER_PSA_LIMITS)]
  
  sim_res = vector("list", nrow(reward_matrix))
  for(i in 1:nrow(reward_matrix)){
    print(paste0("Running for reward row (",i,"): ", paste(reward_matrix[i,], collapse = ' ')))
    sim_res[[i]] = vector("list", length(discount_factors))
    
    for(j in 1:length(discount_factors)){
      DISCOUNT_FACTOR = discount_factors[j]
      print(paste("Running for discount factor:", DISCOUNT_FACTOR))
      sim_res[[i]][[j]] = vector("list", length(max_biopsies_vec))
      
      for(q in 1:length(max_biopsies_vec)){
        max_biopsies = max_biopsies_vec[q]
        print(paste("Running for max_biopsies:", max_biopsies))
        
        sim_res[[i]][[j]][[q]] = vector("list", length(max_depths))
        
        for(max_depth in max_depths){
          
          print(paste("Max depth", max_depth))      
          pat_subset = 1:100
          sim_res[[i]][[j]][[q]][[max_depth + 1]] = testData$testDs.id[pat_subset, c("P_ID", "Age", "progression_time")]
          
          sim_res[[i]][[j]][[q]][[max_depth + 1]][,c('nb', 'offset')] = 
            foreach(pid=sim_res[[i]][[j]][[q]][[max_depth + 1]]$P_ID, 
                    .packages = c("splines", "JMbayes"), .combine = 'rbind') %dopar%{                  
                      #source the common methods for all algorithms
                      source("src/mdp/common/simCommon.R")
                      source("src/mdp/common/prediction.R")
                      source("src/mdp/tree_search/forward_search_no_Y.R")
                      
                      patient_df = testData$testDs[testData$testDs$P_ID==pid,]
                      # if(is.null(patient_df$psa_cat_data)){
                      #   patient_df$psa_cat_data = F
                      # }
                      progression_time = patient_df$progression_time[1]
                      REWARDS = reward_matrix[i,]
                      
                      set.seed(patient_df$P_ID[1])
                      
                      print(paste("Patient ID:", patient_df$P_ID[1], "Progression time:", progression_time))
                      nb = 0
                      delay = -Inf
                      latest_biopsy_time = 0
                      decision_epoch = 1
                      while(decision_epoch<MAX_FOLLOW_UP_TIME & delay<0){
                        pat_data = patient_df[patient_df$visitTimeYears <= decision_epoch,]
                        
                        DISCOUNT_FACTORS = DISCOUNT_FACTOR^(BIOPSY_TEST_TIMES-decision_epoch)
                        names(DISCOUNT_FACTORS) = BIOPSY_TEST_TIMES
                        
                        act = selectAction(pat_data, current_decision_epoch = decision_epoch,
                                           latest_survival_time = latest_biopsy_time, earliest_failure_time = Inf,
                                           max_decision_epoch = min(MAX_FOLLOW_UP_TIME, decision_epoch + max_depth),
                                           cur_biopsies = nb, max_biopsies = max_biopsies)
                        
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
          
          save(sim_res, file = paste0("Rdata/mdp/schedule_by_time/sim_res_",file_num, ".Rdata"))
        }
      }
    }
  }
  rm(fitted_JM)
}
t2 = Sys.time()

stopCluster(ct)