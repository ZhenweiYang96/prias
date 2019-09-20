args = commandArgs(trailingOnly=TRUE)

library(JMbayes)
library(splines)
seed = as.numeric(args[1])

load(paste0("Rdata/lastpaper/sims/sim_seed_", seed, ".Rdata"))

MAX_FOLLOW_UP = 10

source("src/lastpaper/prediction.R")
source("src/lastpaper/minDistThreshold.R")
source("src/clinical_gap3/scheduleCreator.R")
source("src/lastpaper/goodness_of_fit.R")

timesPerSubject = max(sim_res$testData$testDs$visitNumber)
sim_res$testData$testDs$gleason_sum = NA

print("Starting with min distance schedule")

set.seed(as.numeric(seed))
for(nb1_to_delay in c(seq(0,1,0.1), seq(1,2,0.25))){
  
  print(paste('nb1_to_delay is', nb1_to_delay))
  
  nbCol = paste0("mindist_nb_", nb1_to_delay)
  delayCol = paste0("mindist_delay_", nb1_to_delay)
  
  sim_res$testData$testDs.id[, nbCol] = 0
  sim_res$testData$testDs.id[, delayCol] = NA
  
  for(testId in sim_res$testData$testDs.id$P_ID){
    
    print(paste("Doing for patient:", testId))
    curVisitNr = 5
    patient_data = sim_res$testData$testDs[sim_res$testData$testDs$P_ID == testId,]
    patient_data$gleason_sum[1] = 6
    
    idFilter = sim_res$testData$testDs.id$P_ID == testId
    progression_time = sim_res$testData$testDs.id$progression_time[idFilter]
    
    latest_survival_time = 0
    
    while(curVisitNr <= timesPerSubject){
      cur_visit_time = patient_data$year_visit[curVisitNr]
      if(cur_visit_time - latest_survival_time >= 1){
        decision = minDistScheduleDecision(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], 
                                           cur_visit_time = cur_visit_time,
                                           latest_survival_time = latest_survival_time,
                                           nb1_to_delay=nb1_to_delay,
                                           weight_by_horizon_risk = F)
        
        if(decision==T){
          patient_data$gleason_sum[curVisitNr] = 6
          sim_res$testData$testDs.id[idFilter, nbCol] = sim_res$testData$testDs.id[idFilter, nbCol] + 1
          sim_res$testData$testDs.id[idFilter, delayCol] = cur_visit_time - progression_time
          latest_survival_time = cur_visit_time
          if(sim_res$testData$testDs.id[idFilter, delayCol]>=0){
            break;
          }
        }
      }
      
      curVisitNr = curVisitNr + 1
    }
  }
  
  res = sim_res$testData$testDs.id
  save(res, file = paste0("Rdata/lastpaper/res_seed_", seed, ".Rdata"))
}

print("Starting with threshold based schedule")

#Then we do threshold based
for(threshold in seq(0, 1, 0.05)){
  
  print(paste('threshold is', threshold))
  
  nbCol = paste0("threshold_nb_", threshold)
  delayCol = paste0("threshold_delay_", threshold)
  
  sim_res$testData$testDs.id[, nbCol] = 0
  sim_res$testData$testDs.id[, delayCol] = NA
  
  for(testId in sim_res$testData$testDs.id$P_ID){
    
    print(paste("Doing for patient:", testId))
    curVisitNr = 5
    patient_data = sim_res$testData$testDs[sim_res$testData$testDs$P_ID == testId,]
    patient_data$gleason_sum[1] = 6
    
    idFilter = sim_res$testData$testDs.id$P_ID == testId
    progressed = sim_res$testData$testDs.id$progressed[idFilter]
    progression_time = sim_res$testData$testDs.id$progression_time[idFilter]
    
    latest_survival_time = 0
    
    while(curVisitNr <= timesPerSubject){
      cur_visit_time = patient_data$year_visit[curVisitNr]
      
      if(cur_visit_time - latest_survival_time >= 1){
        cum_risk_testpat = 1-getExpectedFutureOutcomes(sim_res$mvJoint_psa_simDs, 
                                                       patient_data[1:curVisitNr,], 
                                                       latest_survival_time, 
                                                       Inf,
                                                       survival_predict_times = cur_visit_time,
                                                       psaDist = "Tdist", M = 500)$predicted_surv_prob
        
        cum_risk_testpat_cur_visit = mean(cum_risk_testpat)
        decision = cum_risk_testpat_cur_visit >= threshold
        
        if(decision==T){
          patient_data$gleason_sum[curVisitNr] = 6
          sim_res$testData$testDs.id[idFilter, nbCol] = sim_res$testData$testDs.id[idFilter, nbCol] + 1
          sim_res$testData$testDs.id[idFilter, delayCol] = cur_visit_time - progression_time
          latest_survival_time = cur_visit_time
          if(progressed & sim_res$testData$testDs.id[idFilter, delayCol]>=0){
            break;
          }
        }
      }
      
      curVisitNr = curVisitNr + 1
    }
  }
  
  res = sim_res$testData$testDs.id
  save(res, file = paste0("Rdata/lastpaper/res_seed_", seed, ".Rdata"))
}

PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP, 0.5))
DRE_CHECK_UP_TIME = seq(0, MAX_FOLLOW_UP, 0.5)
BIOPSY_TEST_TIMES = PSA_CHECK_UP_TIME

roc_results_cache = vector("list", length(BIOPSY_TEST_TIMES))
names(roc_results_cache) = BIOPSY_TEST_TIMES
for(i in 1:length(BIOPSY_TEST_TIMES)){
  roc_results_cache[[i]] = vector("list", length(BIOPSY_TEST_TIMES))
  names(roc_results_cache[[i]]) = BIOPSY_TEST_TIMES
}

print("Starting with FPR based schedule")

for(fpr in seq(0, 1, by = 0.05)){
  
  print(paste('fpr is', fpr))
  
  nbCol = paste0("fpr_nb_", fpr)
  delayCol = paste0("fpr_delay_", fpr)
  
  sim_res$testData$testDs.id[, nbCol] = 0
  sim_res$testData$testDs.id[, delayCol] = NA
  
  for(testId in sim_res$testData$testDs.id$P_ID){
    
    print(paste("Doing for patient:", testId))
    curVisitNr = 5
    patient_data = sim_res$testData$testDs[sim_res$testData$testDs$P_ID == testId,]
    patient_data$gleason_sum[1] = 6
    
    idFilter = sim_res$testData$testDs.id$P_ID == testId
    progressed = sim_res$testData$testDs.id$progressed[idFilter]
    progression_time = sim_res$testData$testDs.id$progression_time[idFilter]
    
    latest_survival_time = 0
    
    while(curVisitNr <= timesPerSubject){
      cur_visit_time = patient_data$year_visit[curVisitNr]
      
      if(cur_visit_time - latest_survival_time >= 1){
        cumrisk_auc_wp = getGaussianQuadWeightsPoints(c(latest_survival_time, cur_visit_time))
        cum_risk_testpat = 1-getExpectedFutureOutcomes(sim_res$mvJoint_psa_simDs, 
                                                       patient_data[1:curVisitNr,], 
                                                       latest_survival_time, 
                                                       Inf,
                                                       survival_predict_times = c(cumrisk_auc_wp$points, cur_visit_time),
                                                       psaDist = "Tdist", M = 500)$predicted_surv_prob
        
        cum_risk_testpat_cur_visit = mean(cum_risk_testpat[nrow(cum_risk_testpat),])
        
        auc_risk_testpat_cur_visit = mean(apply(cum_risk_testpat[-nrow(cum_risk_testpat),],
                                                MARGIN = 2,
                                                FUN = function(risk_col){
                                                  sum(risk_col * cumrisk_auc_wp$weights)
                                                }))
        
        gof = roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]]
        if(is.null(gof)){
          gof = goodness_of_fit(sim_res$mvJoint_psa_simDs, sim_res$trainingData$trainingDs,
                                T_start = latest_survival_time, T_horiz = cur_visit_time, horizon = MAX_FOLLOW_UP, M = 500)
          roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]] = gof
        }
        
        tpr_auc = gof$roc_results_auc_T_start_T_horiz$tpr[which.min(abs(gof$roc_results_auc_T_start_T_horiz$fpr - fpr))]
        tpr_cumrisk = gof$roc_results$tpr[which.min(abs(gof$roc_results$fpr - fpr))]
        if(tpr_auc > tpr_cumrisk){
          threshold = gof$roc_results_auc_T_start_T_horiz$threshold[which.min(abs(gof$roc_results_auc_T_start_T_horiz$fpr - fpr))]
          decision = auc_risk_testpat_cur_visit >= threshold
        }else{
          threshold = gof$roc_results$threshold[which.min(abs(gof$roc_results$fpr - fpr))]
          decision = cum_risk_testpat_cur_visit >= threshold
        }
        
        if(decision==T){
          patient_data$gleason_sum[curVisitNr] = 6
          sim_res$testData$testDs.id[idFilter, nbCol] = sim_res$testData$testDs.id[idFilter, nbCol] + 1
          sim_res$testData$testDs.id[idFilter, delayCol] = cur_visit_time - progression_time
          latest_survival_time = cur_visit_time
          if(progressed & sim_res$testData$testDs.id[idFilter, delayCol]>=0){
            break;
          }
        }
      }
      
      curVisitNr = curVisitNr + 1
    }
  }
  
  res = sim_res$testData$testDs.id
  save(res, file = paste0("Rdata/lastpaper/res_seed_", seed, ".Rdata"))
}

