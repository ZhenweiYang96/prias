library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction.R")
source("src/mdp/tree_search/forward_search_no_Y.R")

max_cores = 3
N_MCMC_ITER = 0
discount_factors = 1
max_biopsies_vec = c(Inf)
max_depths = c(5)

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)
# reward_matrix = as.matrix(expand.grid(seq(5,100,by=10),
#                                       -1, 1,
#                                       seq(-5,-100,by=-10)))
reward_matrix = as.matrix(expand.grid(c(6/12, 12/12, 18/12, 24/12),
                                      -c(0/12, 3/12, 6/12, 12/12), 
                                      c(6/12, 12/12, 18/12, 24/12),
                                      -1))

colnames(reward_matrix) = reward_names

ct = makeCluster(max_cores, type = "FORK")
registerDoParallel(ct)

nDs = 10
for(file_num in 10){
  
  file = list.files(path='/home/a_tomer/Data/final_res_2nd_paper/', full.names = T)[file_num]
  load(file)
  
  print(paste("Loaded file:", file_num))
  
  fitted_JM = jointModelData$mvJoint_dre_psa_simDs
  jointModelData$mvJoint_dre_psa_simDs = NULL
  
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
          pidl = 1:100
          sim_res[[i]][[j]][[q]][[max_depth + 1]] = jointModelData$testData$testDs.id[pidl, c("P_ID", "Age", "progression_time")]
          
          sim_res[[i]][[j]][[q]][[max_depth + 1]][,c('nb', 'offset')] = 
            foreach(pid=sim_res[[i]][[j]][[q]][[max_depth + 1]]$P_ID, 
                    .packages = c("splines", "JMbayes"), .combine = 'rbind') %dopar%{                  
                      patient_df = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID==pid,]
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
                                           G6_probs = NULL, latest_survival_time = latest_biopsy_time, earliest_failure_time = Inf,
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

stopCluster(ct)
