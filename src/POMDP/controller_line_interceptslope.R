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
N_DESPOT_SCENARIOS = 500
discount_factors = c(1, 0.5, 0)
max_depth = 4
max_biopsies = Inf

file = list.files(path='/home/a_tomer/Data/final_res_2nd_paper/', full.names = T)[1]
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

#These two are intercepts and slope of wait (angles in vector of reward_fw)
angles = c(1, 5, 15, 45, 60, 80)
reward_tw = 0
rewards_fw = -tan(pi/180 * angles) * 1

thresholds = c(0.05, 0.10, 0.15, 0.30, 0.50)

ct = makeCluster(max_cores)
registerDoParallel(ct)

sim_res = vector("list", length(thresholds))
names(sim_res)=thresholds
#for(t in 1:length(thresholds)){
for(t in c(2,4,1,5,3)){
  threshold = thresholds[t]
  print(paste("Running for threshold number",t, ":", threshold))
  
  sim_res[[t]] = vector("list", length(discount_factors))
  names(sim_res[[t]])=discount_factors
  
  for(df in 1:length(discount_factors)){
    DISCOUNT_FACTOR = discount_factors[df]
    print(paste("Running for discount_factor number",df, ":", DISCOUNT_FACTOR))
    
    sim_res[[t]][[df]] = vector("list", length(rewards_fw))
    names(sim_res[[t]][[df]])=rewards_fw
    
    for(rfw in 1:length(rewards_fw)){
      reward_fw = rewards_fw[rfw]
      print(paste("Running for false wait reward number", rfw, ":", reward_fw))
      
      dist = threshold * reward_fw
      rewards_fb = seq(dist/2, dist - threshold / tan(pi/180 * 25), length.out=5)
      
      slopes = (dist - rewards_fb) / threshold
      rewards_tb = rewards_fb + slopes * 1
      
      sim_res[[t]][[df]][[rfw]] = vector("list", length(rewards_fb))
      
      for(rtbfb in 1:length(rewards_fb)){
        reward_fb = rewards_fb[rtbfb]
        reward_tb = rewards_tb[rtbfb]
        
        print(paste("Running for truebiopsy/falsebiopsy number", rtbfb)) 
        print(paste("Slope of biopsy line: ", slopes[rtbfb]))
        print(paste("Reward false biopsy: ", reward_fb)) 
        print(paste("Reward true biopsy: ", reward_tb)) 
        
        reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)
        REWARDS = c(reward_tb, reward_fb, reward_tw, reward_fw)
        names(REWARDS) = reward_names
        
        pat_subset = 1:100
        sim_res[[t]][[df]][[rfw]][[rtbfb]] = testData$testDs.id[pat_subset, c("P_ID", "Age", "progression_time")]
        
        sim_res[[t]][[df]][[rfw]][[rtbfb]][,c('nb', 'offset')] = 
          foreach(pid=sim_res[[t]][[df]][[rfw]][[rtbfb]]$P_ID, 
                  .packages = c("splines", "JMbayes"), 
                  .combine = 'rbind') %dopar%{
                    
                    source("src/mdp/common/simCommon.R")
                    source("src/mdp/common/prediction.R")
                    source("src/mdp/tree_search/forward_search_no_Y.R")
                    
                    patient_df = testData$testDs[testData$testDs$P_ID==pid,]
                    # if(is.null(patient_df$psa_cat_data)){
                    #   patient_df$psa_cat_data = F
                    # }
                    progression_time = patient_df$progression_time[1]
                    
                    set.seed(pid)
                    
                    print(paste("Patient ID:", pid, "Progression time:", progression_time))
                    nb = 0
                    delay = -Inf
                    latest_biopsy_time = 0
                    decision_epoch = 1
                    
                    while(decision_epoch<MAX_FOLLOW_UP_TIME & delay<0){
                      pat_data = patient_df[patient_df$visitTimeYears <= decision_epoch,]
                      
                      DISCOUNT_FACTORS = DISCOUNT_FACTOR^(BIOPSY_TEST_TIMES-decision_epoch)
                      names(DISCOUNT_FACTORS) = BIOPSY_TEST_TIMES
                      
                      if(DISCOUNT_FACTOR!=0){
                        max_decision_epoch = min(MAX_FOLLOW_UP_TIME, decision_epoch + max_depth)
                      }else{
                        max_decision_epoch = decision_epoch
                      }
                      
                      act = selectAction(pat_data, current_decision_epoch = decision_epoch,
                                         latest_survival_time = latest_biopsy_time, earliest_failure_time = Inf,
                                         max_decision_epoch = max_decision_epoch,
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
        res = sim_res[[t]]
        save(res, file = paste0("Rdata/mdp/schedule_by_reward/sim_res",
                                threshold, ".Rdata"))
      }
    }
  }
}
stopCluster(ct)
