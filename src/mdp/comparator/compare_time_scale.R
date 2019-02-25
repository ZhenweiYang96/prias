setwd("/home/atomer/src/prias/")

library(JMbayes)
library(survival)
library(splines)
library(ggplot2)
library(doParallel)

#Load a simulation
load("Rdata/mdp/jointModelData_seed_101_simNr_1_normal.Rdata")

fitted_JM = jointModelData$mvJoint_dre_psa_simDs
jointModelData$mvJoint_dre_psa_simDs = NULL

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
#Source the method you want to use
source("src/mdp/tree_search/forward_search_no_Y.R")

max_cores = 10

decision_epochs = 1:10
max_decision_epochs = sapply(decision_epochs, FUN = function(de){
  seq(de, min(de + 4, 10), by = 1)
}, simplify = F)
names(max_decision_epochs) = as.character(decision_epochs)

time_to_biopsy_scales = unique(c(10/0.5, 7/0.7, 5/1.1, 4/1.8, 4/0.7, 3/0.8))
names(time_to_biopsy_scales) = c("Annual", "1pt5years", "Biennial", "Triennial",
                                "PRIAS_risk10pc", "risk15pc")

ct = makeCluster(max_cores)
registerDoParallel(ct)
sim_res = vector("list", length(time_to_biopsy_scales))
names(sim_res) = names(time_to_biopsy_scales)
for(i in 1:length(time_to_biopsy_scales)){
  sim_res[[i]] = vector("list", length(decision_epochs))
  names(sim_res[[i]]) = as.character(decision_epochs)
  
  for(decision_epoch in decision_epochs){
    decision_epoch_str = as.character(decision_epoch)
    
    dataset = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears<=decision_epoch,]
    dataset = split(dataset, dataset$P_ID)
    
    sim_res[[i]][[decision_epoch_str]] = vector("list", length(max_decision_epochs[[decision_epoch_str]]))
    names(sim_res[[i]][[decision_epoch_str]]) = as.character(max_decision_epochs[[decision_epoch_str]])
    
    for(max_decision_epoch in max_decision_epochs[[decision_epoch_str]]){
      max_decision_epoch_str = as.character(max_decision_epoch)
      
      sim_res[[i]][[decision_epoch_str]][[max_decision_epoch_str]] = 
        foreach(j=1:length(dataset), 
                .export=c("fitted_JM"),
                .packages = c("splines", "JMbayes")) %dopar%{
                  source("src/mdp/common/simCommon.R")
                  source("src/mdp/tree_search/forward_search_no_Y.R")
                  
                  patient_df = dataset[[j]]
                  time_to_biopsy_scale = time_to_biopsy_scales[i]
                  #REWARDS = reward_matrix[i,]
                  
                  DISCOUNT_FACTORS = DISCOUNT_FACTOR^((1:length(PSA_CHECK_UP_TIME))-nrow(patient_df))
                  names(DISCOUNT_FACTORS) = PSA_CHECK_UP_TIME
                  
                  set.seed(patient_df$P_ID[1])
                  selectAction(patient_df, current_decision_epoch = decision_epoch,
                               latest_survival_time = 0, earliest_failure_time = Inf,
                               max_decision_epoch = max_decision_epoch)            
                }
      
      save(sim_res, file = paste0("Rdata/mdp/reward_by_time/sim_res_", 
                                  names(time_to_biopsy_scales)[i], ".Rdata"))
    }
  }
}
stopCluster(ct)
