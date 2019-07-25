rm(list = ls())
load("Rdata/lastpaper/sim_seed_2019.Rdata")
source("src/lastpaper/prediction.R")
source("src/lastpaper/goodness_of_fit.R")

timesPerSubject = max(sim_res$testData$testDs$visitNumber)

sim_res$testData$testDs$gleason_sum = NA

PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP, 0.5))
DRE_CHECK_UP_TIME = seq(0, MAX_FOLLOW_UP, 0.5)
BIOPSY_TEST_TIMES = DRE_CHECK_UP_TIME

roc_results_cache = vector("list", length(BIOPSY_TEST_TIMES))
names(roc_results_cache) = BIOPSY_TEST_TIMES
lapply(1:length(BIOPSY_TEST_TIMES), function(i){
  roc_results_cache[[i]] <<- vector("list", length(BIOPSY_TEST_TIMES))
  names(roc_results_cache[[i]]) <<- BIOPSY_TEST_TIMES
})

##################################################
#FPR based schedule. both auc_risk threshold and cumrisk
##################################################
for(fpr in seq(0, 0.35, by = 0.05)){
  
  nbCol = paste0("nb_", fpr)
  delayCol = paste0("delay_", fpr)
  
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
        set.seed(2019)
        cumrisk_auc_wp = getGaussianQuadWeightsPoints(c(latest_survival_time, cur_visit_time))
        cum_risk_testpat = 1-getExpectedFutureOutcomes(sim_res$mvJoint_psa_simDs, 
                                                       patient_data[1:curVisitNr,], 
                                                       latest_survival_time, 
                                                       Inf,
                                                       survival_predict_times = c(cumrisk_auc_wp$points, cur_visit_time),
                                                       psaDist = "Tdist", M = 800)$predicted_surv_prob
        
        cum_risk_testpat_cur_visit = mean(cum_risk_testpat[nrow(cum_risk_testpat),])
        
        auc_risk_testpat_cur_visit = mean(apply(cum_risk_testpat[-nrow(cum_risk_testpat),],
                                                MARGIN = 2,
                                                FUN = function(risk_col){
                                                  sum(risk_col * cumrisk_auc_wp$weights)
                                                }))
        
        gof = roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]]
        if(is.null(gof)){
          gof = goodness_of_fit(sim_res$mvJoint_psa_simDs, sim_res$trainingData$trainingDs,
                                T_start = latest_survival_time, T_horiz = cur_visit_time, horizon = MAX_FOLLOW_UP, M = 800)
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
  
  testDs.id = sim_res$testData$testDs.id
  save(testDs.id, file=paste0("Rdata/lastpaper/fpr/sim_aucrisk_cumrisk_fpr_", fpr,".Rdata"))
}

##################################################
#FPR based schedule. Using only auc_risk threshold
##################################################
for(fpr in seq(0, 0.35, by = 0.05)){
  
  nbCol = paste0("nb_", fpr)
  delayCol = paste0("delay_", fpr)
  
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
        set.seed(2019)
        cumrisk_auc_wp = getGaussianQuadWeightsPoints(c(latest_survival_time, cur_visit_time))
        cum_risk_testpat = 1-getExpectedFutureOutcomes(sim_res$mvJoint_psa_simDs, 
                                                       patient_data[1:curVisitNr,], 
                                                       latest_survival_time, 
                                                       Inf,
                                                       survival_predict_times = cumrisk_auc_wp$points,
                                                       psaDist = "Tdist", M = 800)$predicted_surv_prob
        
        auc_risk_testpat_cur_visit = mean(apply(cum_risk_testpat, MARGIN = 2,
                                                FUN = function(risk_col){
                                                  sum(risk_col * cumrisk_auc_wp$weights)
                                                }))
        
        gof = roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]]
        if(is.null(gof)){
          gof = goodness_of_fit(sim_res$mvJoint_psa_simDs, sim_res$trainingData$trainingDs,
                                T_start = latest_survival_time, T_horiz = cur_visit_time, horizon = MAX_FOLLOW_UP, M = 800)
          roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]] = gof
        }
        
        threshold = gof$roc_results_auc_T_start_T_horiz$threshold[which.min(abs(gof$roc_results_auc_T_start_T_horiz$fpr - fpr))]
        decision = auc_risk_testpat_cur_visit >= threshold
        
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
  
  testDs.id = sim_res$testData$testDs.id
  save(testDs.id, file=paste0("Rdata/lastpaper/fpr/sim_aucrisk_fpr_", fpr,".Rdata"))
}


##################################################
#FPR based schedule. Using only cum risk threshold
##################################################
for(fpr in seq(0, 0.35, by = 0.05)){
  
  nbCol = paste0("nb_", fpr)
  delayCol = paste0("delay_", fpr)
  
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
        set.seed(2019)
        cum_risk_testpat = 1-getExpectedFutureOutcomes(sim_res$mvJoint_psa_simDs, 
                                                       patient_data[1:curVisitNr,], 
                                                       latest_survival_time, 
                                                       Inf,
                                                       survival_predict_times = cur_visit_time,
                                                       psaDist = "Tdist", M = 800)$predicted_surv_prob
        
        cum_risk_testpat_cur_visit = mean(cum_risk_testpat)
        
        gof = roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]]
        if(is.null(gof)){
          gof = goodness_of_fit(sim_res$mvJoint_psa_simDs, sim_res$trainingData$trainingDs,
                                T_start = latest_survival_time, T_horiz = cur_visit_time, horizon = MAX_FOLLOW_UP, M = 800)
          roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]] = gof
        }
        
        threshold = gof$roc_results$threshold[which.min(abs(gof$roc_results$fpr - fpr))]
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
  
  testDs.id = sim_res$testData$testDs.id
  save(testDs.id, file=paste0("Rdata/lastpaper/fpr/sim_cumrisk_fpr_", fpr,".Rdata"))
}