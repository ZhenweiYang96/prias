library(JMbayes)
library(survival)
library(splines)
library(ggplot2)
library(doParallel)

#Load a simulation
load("/home/a_tomer/Data/final_res_2nd_paper/jointModelData_seed_101_simNr_1_normal.Rdata")
#load("C:/Users/838035/jointModelData_seed_101_simNr_1_normal.Rdata")

fitted_JM = jointModelData$mvJoint_dre_psa_simDs
jointModelData$mvJoint_dre_psa_simDs = NULL

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
#Source the method you want to use
source("src/mdp/tree_search/forward_search_no_Y.R")

testPatientIds = jointModelData$testData$testDs.id$P_ID

max_cores = 3

#for(threshold in c(0.15, 0.125, 0.1, 0.075, 0.05)){
#for(threshold in c(0.10)){
  ct = makeCluster(max_cores)
  registerDoParallel(ct)
  t1 = Sys.time()
  test_times = 1
  
  dataset.id = comparison_res = vector("list", length(test_times))
  names(dataset.id) = names(comparison_res) = test_times
  for(test_time in test_times){
    test_time_str = as.character(test_time)
    
    dataset = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears<=test_time,]
    dataset.id[[test_time_str]] = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears==test_time,]
    
    dataset = dataset[dataset$P_ID %in% testPatientIds,]
    dataset = split(dataset, dataset$P_ID)
    
    max_decision_epochs = c(test_time, 1.25, 1.5, 1.75, 2, 2.5,3,3.5,4,4.5,5)
    #max_decision_epochs = c(test_time, 1.25, 1.75, 2, 3,4,5)
    comparison_res[[test_time_str]] = vector("list", length(max_decision_epochs))
    names(comparison_res[[test_time_str]]) = max_decision_epochs
    
    for(max_decision_epoch in max_decision_epochs){
      max_decision_epoch_str = as.character(max_decision_epoch)
      dataset.id[[test_time_str]][,max_decision_epoch_str] = NA
      comparison_res[[test_time_str]][[max_decision_epoch_str]] = foreach(i=1:length(dataset), .export=c("fitted_JM"),
                                                                          .packages = c("splines", "JMbayes")) %dopar%{
                                                                            patient_df = dataset[[i]]
                                                                            
                                                                            # set.seed(patient_df$P_ID[1] + 900)
                                                                            # slope_b = runif(1, 0, 10)
                                                                            # intercept_b = runif(1, -10, 10)
                                                                            # 
                                                                            # intercept_w = runif(1, intercept_b + 1, intercept_b + 10)
                                                                            # slope_w = (intercept_b + slope_b * threshold - intercept_w)/threshold
                                                                            
                                                                            #source the common methods for all algorithms
                                                                            source("src/mdp/common/simCommon.R")
                                                                            source("src/mdp/tree_search/forward_search_no_Y.R")
                                                                            
                                                                            DISCOUNT_FACTORS = DISCOUNT_FACTOR^((1:length(PSA_CHECK_UP_TIME))-nrow(patient_df))
                                                                            names(DISCOUNT_FACTORS) = PSA_CHECK_UP_TIME
                                                                            
                                                                            set.seed(patient_df$P_ID[1])
                                                                            selectAction(patient_df, current_decision_epoch = test_time,
                                                                                         latest_survival_time = 0, earliest_failure_time = Inf,
                                                                                         max_decision_epoch = max_decision_epoch)            
                                                                          }
      
      decisions = sapply(comparison_res[[test_time_str]][[max_decision_epoch_str]],FUN = function(x){x$optimal_action})
      dataset.id[[test_time_str]][,max_decision_epoch_str] = decisions
    }
  }
  t2 = Sys.time()
  t2 - t1
  stopCluster(ct)
  
  #save(dataset.id, file = paste0("Rdata/mdp/DF_pt9/result_", threshold, ".Rdata"))
#}
