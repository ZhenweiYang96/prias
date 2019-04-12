library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction.R")
source("src/mdp/tree_search/forward_search_no_Y.R")

max_cores = 3
N_MCMC_ITER = 250
N_DESPOT_SCENARIOS = 500
discount_factors = c(0, 1, 0.5)
max_depth = 3

file = list.files(path='C:/Users/838035/prias/Rdata/mdp/final_res/', full.names = T)[1]
load(file)
print(paste("Loaded file:", 1))

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



value_functions = c(10)
thresholds = c(0.15, 0.10)

#slopes = c(0.414, 1, 2.414, 0, -0.414, -1, -2.414)
slopes = c(100, 1, 0.2, 0, -1, -0.2, -100)

ct = makeCluster(max_cores)
registerDoParallel(ct)

sim_res = vector("list", length(thresholds))
names(sim_res)=thresholds
for(t in 1:length(thresholds)){
  threshold = thresholds[t]
  print(paste("Running for threshold number",t, ":", threshold))
  
  sim_res[[t]] = vector("list", length(discount_factors))
  names(sim_res[[t]])=discount_factors
  
  for(df in 1:length(discount_factors)){
    DISCOUNT_FACTOR = discount_factors[df]
    print(paste("Running for discount_factor number",df, ":", DISCOUNT_FACTOR))
    
    sim_res[[t]][[df]] = vector("list", length(value_functions))
    names(sim_res[[t]][[df]])=value_functions
    
    for(v in 1:length(value_functions)){
      value = value_functions[v]
      print(paste("Running for value number",v, ":", value))
      
      sim_res[[t]][[df]][[v]] = vector("list", length(slopes))
      names(sim_res[[t]][[df]][[v]])=slopes
      
      for(bs in 1:length(slopes)){
        biopsy_slope = slopes[bs]
        print(paste("Running for biopsy slope number",bs, ":", biopsy_slope))
        
        biopsy_intercept = value - biopsy_slope * threshold
        print(paste("Biopsy intercept: ", biopsy_intercept))
        
        sim_res[[t]][[df]][[v]][[bs]] = vector("list", length(slopes))
        names(sim_res[[t]][[df]][[v]][[bs]])=slopes
        
        for(ws in 1:length(slopes)){
          wait_slope = slopes[ws]
          
          if(wait_slope != biopsy_slope){
            print(paste("Running for wait slope number",ws, ":", wait_slope))
            wait_intercept = value - wait_slope*threshold
            print(paste("Wait intercept: ", wait_intercept))
            
            REWARDS = thresholdToReward(NA, biopsy_intercept, biopsy_slope, wait_intercept, wait_slope)
            print(REWARDS)
            
            pat_subset = 1:100
            sim_res[[t]][[df]][[v]][[bs]][[ws]] = testData$testDs.id[pat_subset, c("P_ID", "Age", "progression_time")]
            
            sim_res[[t]][[df]][[v]][[bs]][[ws]][,c('optimal_action', 
                                                   'optimal_action_chain',
                                                   'optimal_reward')] = 
              foreach(pid=sim_res[[t]][[df]][[v]][[bs]][[ws]]$P_ID, 
                      .packages = c("splines", "JMbayes"), 
                      .combine = 'rbind') %dopar%{
                        source("src/mdp/common/simCommon.R")
                        source("src/mdp/common/prediction.R")
                        source("src/mdp/tree_search/forward_search_no_Y.R")
                        
                        patient_df = testData$testDs[testData$testDs$P_ID==pid,]
                        if(is.null(patient_df$psa_cat_data)){
                          patient_df$psa_cat_data = F
                        }
                        progression_time = patient_df$progression_time[1]
                        
                        set.seed(pid)
                        
                        print(paste("Patient ID:", pid, "Progression time:", progression_time))
                        nb = 0
                        delay = -Inf
                        latest_biopsy_time = 0
                        decision_epoch = 2
                        pat_data = patient_df[patient_df$visitTimeYears <= decision_epoch,]
                        
                        DISCOUNT_FACTORS = DISCOUNT_FACTOR^(BIOPSY_TEST_TIMES-decision_epoch)
                        names(DISCOUNT_FACTORS) = BIOPSY_TEST_TIMES
                        
                        max_decision_epoch = if(DISCOUNT_FACTOR==0){
                          decision_epoch
                        }else{
                          decision_epoch + max_depth
                        }
                        
                        act = selectAction(pat_data, current_decision_epoch = decision_epoch,
                                           latest_survival_time = latest_biopsy_time, earliest_failure_time = Inf,
                                           max_decision_epoch = max_decision_epoch,
                                           cur_biopsies = nb, max_biopsies = Inf)
                        
                        print(act$G6_probs[1])
                        return(c('optimal_action'=act$optimal_action, 
                                 'optimal_action_chain'=act$optimal_action_chain,
                                 'optimal_reward'=act$optimal_reward))
                      }
            save(sim_res, file = paste0("Rdata/mdp/decision_for_lines/sim_res_",
                                        file_num, ".Rdata"))
          }
        }
      }
    }
  }
}

stopCluster(ct)