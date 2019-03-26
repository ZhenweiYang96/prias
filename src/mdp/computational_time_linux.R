setwd("/home/a_tomer/Google Drive/PhD/src/prias")

library(pryr)
library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction.R")
source("src/mdp/tree_search/forward_search_no_Y.R")

N_MCMC_ITER = 500

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)
REWARDS = c(100, -1, 1, -100)
names(REWARDS) = reward_names

m1.load = mem_used()

t1.load = Sys.time()
file = list.files(path='/home/a_tomer/Data/final_res_2nd_paper', full.names = T)[9]
load(file)
t2.load = Sys.time()

m2.load = mem_used()

print(t2.load - t1.load)

fitted_JM = jointModelData$mvJoint_dre_psa_simDs
jointModelData$mvJoint_dre_psa_simDs = NULL

DISCOUNT_FACTOR = 1
max_biopsies = Inf
max_depth = 5        

sim_res = jointModelData$testData$testDs.id[1:50, c("P_ID", "Age", "progression_time")]

m1Common.sim = mem_used()
t1Common.sim = Sys.time()

ct = makeCluster(3, type='FORK')
registerDoParallel(ct)

sim_res[,c('t1', 't2','m1', 'm2', 'nb')] = foreach(pid=sim_res$P_ID, 
                                      .packages = c("splines", "JMbayes", "pryr"), .combine = 'rbind') %dopar%{                  
                                        
                                        source("src/mdp/common/simCommon.R")
                                        source("src/mdp/common/prediction.R")
                                        source("src/mdp/tree_search/forward_search_no_Y.R")
                                        
                                        m1.sim = mem_used()
                                        t1.sim = Sys.time()
                                        patient_df = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID==pid,]
                                        progression_time = patient_df$progression_time[1]
                                        
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
                                                             cur_biopsies = 0, max_biopsies = max_biopsies)
                                          
                                          if(act$optimal_action==BIOPSY){
                                            latest_biopsy_time = decision_epoch
                                            delay = decision_epoch - progression_time
                                            decision_epoch = BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES >= (decision_epoch + MIN_BIOPSY_GAP)][1]
                                            nb = nb + 1
                                          }else{
                                            decision_epoch = getNextDecisionEpoch(decision_epoch)
                                          }
                                        }
                                        
                                        t2.sim = Sys.time()
                                        m2.sim = mem_used()
                                        
                                        return(c('t1'=t1.sim, 't2'=t2.sim,
                                                 'm1'=m1.sim, 'm2'=m2.sim, 'nb'=nb))
                                      }

t2Common.sim = Sys.time()
m2Common.sim = mem_used()

stopCluster(ct)
