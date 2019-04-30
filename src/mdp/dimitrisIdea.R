library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction.R")

max_cores = 3
N_MCMC_ITER = 250
N_DESPOT_SCENARIOS = 500
discount_factors = c(1, 0.5)
max_depth = 3

file = list.files(path='C://Users//838035//prias//Rdata//mdp//final_res', full.names = T)[1]
load(file)
print(paste("Loaded file:", 1))

fitted_JM = jointModelData$mvJoint_dre_psa_simDs
testData = jointModelData$testData
jointModelData$mvJoint_dre_psa_simDs = NULL

thresholds = c(0.15, 0.10)

ct = makeCluster(max_cores)
registerDoParallel(ct)

sim_res = vector("list", length(thresholds))
names(sim_res)=thresholds
for(t in 1:length(thresholds)){
  threshold = thresholds[t]
  print(paste("Running for threshold number",t, ":", threshold))
  
  sim_res[[t]] = testData$testDs.id[, c("P_ID", "Age", "progression_time")]
  
  biopsyres_colnames = c('B1_prob', 'B1_time', 'B1_delay',
               'B2_prob', 'B2_time', 'B2_delay',
               'B3_prob', 'B3_time', 'B3_delay')
  sim_res[[t]][,biopsyres_colnames] = 
    foreach(pid=sim_res[[t]]$P_ID, 
            .packages = c("splines", "JMbayes"), 
            .combine = 'rbind') %dopar%{
              
              source("src/mdp/common/simCommon.R")
              source("src/mdp/common/prediction.R")
              
              patient_df = testData$testDs[testData$testDs$P_ID==pid,]
              progression_time = patient_df$progression_time[1]
              
              set.seed(pid)
              
              print(paste("Patient ID:", pid, "Progression time:", progression_time))
              
              nb = 0
              delay = -Inf
              latest_biopsy_time = 0
              decision_epoch = 1
              
              biopsyres = rep(NA,9)
              names(biopsyres) = biopsyres_colnames
              
              pat_data = patient_df[patient_df$visitTimeYears <= decision_epoch,]
              
              curCumRisk = NULL
              while(decision_epoch < MAX_FOLLOW_UP_TIME & nb < 3){
                if((decision_epoch - latest_biopsy_time) >= MIN_BIOPSY_GAP){
                  
                  if(is.null(curCumRisk)){
                    curCumRisk = 1 - rowMeans(getExpectedFutureOutcomes(fitted_JM, pat_data, latest_biopsy_time, Inf, 
                                                             BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES>=decision_epoch], M = N_MCMC_ITER)$predicted_surv_prob)
                  }
                  
                  if(curCumRisk[1] > threshold){
                    wtPoints = getGaussianQuadWeightsPoints(c(latest_biopsy_time, decision_epoch))
                    points = sort(wtPoints$points, decreasing = F)
                    wt = wtPoints$weights[order(wtPoints$points, decreasing = F)]
                    
                    survProbMatrix = getExpectedFutureOutcomes(fitted_JM, pat_data, latest_survival_time = latest_biopsy_time, 
                                                               earliest_failure_time = decision_epoch, 
                                                               points, M = N_MCMC_ITER)$predicted_surv_prob
                    latest_biopsy_time = decision_epoch
                    expectedDelay = latest_biopsy_time - sum(wt*rowMeans(survProbMatrix))
                    
                    biopsyres[(nb*3 + 1):(nb*3 + 3)] = c(curCumRisk[1], latest_biopsy_time, expectedDelay)
                    nb = nb + 1
                    
                    curCumRisk = NULL
                  }else{
                    curCumRisk = curCumRisk[-1]
                  }
                }
                
                decision_epoch = getNextDecisionEpoch(decision_epoch)
              }
              
              return(biopsyres)
            }
}

stopCluster(ct)