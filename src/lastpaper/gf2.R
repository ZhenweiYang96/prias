goodness_of_fit <- function (object, newdata, T_start, T_horiz, horizon = 10, M = 500) {
  
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
  #DRE check up time years
  DRE_CHECK_UP_TIME = seq(0, horizon, 0.5)
  BIOPSY_TEST_TIMES = DRE_CHECK_UP_TIME
  
  #First make sure that longitudinal data is available only upto T_horiz
  newdata = newdata[newdata$year_visit <= T_horiz,]
  
  #Now select patients who have a progression time more than T_start
  #This automatically means that they had no event until T_start
  newdata = newdata[newdata$right_cens_time > T_start,]
  newdata.id = newdata[!duplicated(newdata$P_ID),]
  
  ########################
  #Make sure there are patients to calculate TPR etc
  ########################
  
  #Now lets start with real patient status
  newdata.id$real_period_status = rep(NA, nrow(newdata.id))
  newdata.id$real_period_status[newdata.id$reclassification==T & newdata.id$right_cens_time <= T_horiz] = 1
  newdata.id$real_period_status[newdata.id$reclassification==T & newdata.id$right_cens_time > T_horiz] = 0
  newdata.id$real_period_status[newdata.id$reclassification==F & newdata.id$right_cens_time >= T_horiz] = 0
  
  #Now for patients who had reclassification = F and right_cens_time <= T_horiz
  # Their latest biopsy time is same right_cens_time
  # we need to calculate their cum risk using their real data
  subset_patients = newdata[newdata$reclassification==F & newdata$right_cens_time < T_horiz,]
  
  if(nrow(subset_patients) > 0){
    subset_patients_cum_risk_Thoriz = by(INDICES = subset_patients$P_ID, data = subset_patients, FUN = function(pat_data){
      rowMeans(1 - getExpectedFutureOutcomes(object, pat_data, 
                                             latest_survival_time = pat_data$right_cens_time[1],
                                             earliest_failure_time = Inf,
                                             survival_predict_times = T_horiz,
                                             psaDist = "Tdist", M=M)$predicted_surv_prob)
    })
    
    newdata.id$real_period_status[newdata.id$reclassification==F & 
                                    newdata.id$right_cens_time < T_horiz] = as.numeric(subset_patients_cum_risk_Thoriz)
  }
  
  #Now for each of the patients obtain their cumulative risk at 
  #various follow-ups under the condition
  #that they didnt have an event until T_start
  cumrisk_auc_wp = getGaussianQuadWeightsPoints(c(T_start, T_horiz))
  cumrisk_auc_wp_2 = getGaussianQuadWeightsPoints(c(T_horiz, T_horiz + 1))
  cum_risk_T_start_onwards = by(INDICES = newdata$P_ID, data = newdata, FUN = function(pat_data){
    1 - getExpectedFutureOutcomes(object, pat_data, 
                                  latest_survival_time = T_start,
                                  earliest_failure_time = Inf,
                                  survival_predict_times = c(cumrisk_auc_wp$points, T_horiz, cumrisk_auc_wp_2$points),
                                  psaDist = "Tdist", M = M)$predicted_surv_prob
  })
  
  newdata.id$cum_risk_T_start_T_horiz = as.numeric(sapply(cum_risk_T_start_onwards, function(x){
    mean(x[nrow(x),])
  }))
  
  newdata.id$auc_T_start_T_horiz = as.numeric(sapply(cum_risk_T_start_onwards, function(x){
    aucs = apply(x[1:15,], MARGIN = 2, FUN = function(risk_col){
      sum(risk_col * cumrisk_auc_wp$weights)
    })
    
    return(mean(aucs))
  }))
  
  newdata.id$auc_T_start_T_horiz1 = as.numeric(sapply(cum_risk_T_start_onwards, function(x){
    aucs1 = apply(x[1:15,], MARGIN = 2, FUN = function(risk_col){
      sum(risk_col * cumrisk_auc_wp$weights)
    })
    
    aucs2 = apply(x[17:31,], MARGIN = 2, FUN = function(risk_col){
      sum(risk_col * cumrisk_auc_wp_2$weights)
    })
    
    return(mean(aucs1 + aucs2))
  }))
  
  ###############
  # First we get normal ROC results
  ################
  minthres = max(0, min(newdata.id$cum_risk_T_start_T_horiz) - 0.001)
  maxthres = min(1, max(newdata.id$cum_risk_T_start_T_horiz) + 0.001)
  thresholds_cumrisk = seq(minthres, maxthres, length.out = 500)
  
  roc_results = t(sapply(thresholds_cumrisk, FUN = function(threshold){
    predicted_cancer = newdata.id$cum_risk_T_start_T_horiz >= threshold
    real_cancer = newdata.id$real_period_status
    
    c("threshold" = threshold,
      "nTP" = sum(real_cancer * predicted_cancer),
      "nFN" = sum(real_cancer * (1-predicted_cancer)),
      "nFP" = sum((1-real_cancer) * predicted_cancer),
      "nTN" = sum((1-real_cancer) * (1-predicted_cancer)))
  }))
  
  roc_results = data.frame(roc_results)
  roc_results$tpr = roc_results$nTP / (roc_results$nTP + roc_results$nFN)
  roc_results$fpr = roc_results$nFP / (roc_results$nTN + roc_results$nFP)
  
  ###############
  # Then we get auc_T_start_T_horiz ROC results
  ################
  minthres = max(0, min(newdata.id$auc_T_start_T_horiz) - 0.001)
  maxthres = min(horizon, max(newdata.id$auc_T_start_T_horiz) + 0.001)
  thresholds_aucrisk = seq(minthres, maxthres, length.out = 500)
  
  roc_results_auc_T_start_T_horiz = t(sapply(thresholds_aucrisk, FUN = function(threshold){
    predicted_cancer = newdata.id$auc_T_start_T_horiz >= threshold
    real_cancer = newdata.id$real_period_status
    
    c("threshold" = threshold,
      "nTP" = sum(real_cancer * predicted_cancer),
      "nFN" = sum(real_cancer * (1-predicted_cancer)),
      "nFP" = sum((1-real_cancer) * predicted_cancer),
      "nTN" = sum((1-real_cancer) * (1-predicted_cancer)))
  }))
  
  roc_results_auc_T_start_T_horiz = data.frame(roc_results_auc_T_start_T_horiz)
  roc_results_auc_T_start_T_horiz$tpr = roc_results_auc_T_start_T_horiz$nTP / (roc_results_auc_T_start_T_horiz$nTP + roc_results_auc_T_start_T_horiz$nFN)
  roc_results_auc_T_start_T_horiz$fpr = roc_results_auc_T_start_T_horiz$nFP / (roc_results_auc_T_start_T_horiz$nTN + roc_results_auc_T_start_T_horiz$nFP)
  
  ###############
  # Then we get auc_T_start_T_horiz1 ROC results
  ################
  minthres = max(0, min(newdata.id$auc_T_start_T_horiz1) - 0.001)
  maxthres = min(horizon, max(newdata.id$auc_T_start_T_horiz1) + 0.001)
  thresholds_aucrisk = seq(minthres, maxthres, length.out = 500)
  
  roc_results_auc_T_start_T_horiz1 = t(sapply(thresholds_aucrisk, FUN = function(threshold){
    predicted_cancer = newdata.id$auc_T_start_T_horiz1 >= threshold
    real_cancer = newdata.id$real_period_status
    
    c("threshold" = threshold,
      "nTP" = sum(real_cancer * predicted_cancer),
      "nFN" = sum(real_cancer * (1-predicted_cancer)),
      "nFP" = sum((1-real_cancer) * predicted_cancer),
      "nTN" = sum((1-real_cancer) * (1-predicted_cancer)))
  }))
  
  roc_results_auc_T_start_T_horiz1 = data.frame(roc_results_auc_T_start_T_horiz1)
  roc_results_auc_T_start_T_horiz1$tpr = roc_results_auc_T_start_T_horiz1$nTP / (roc_results_auc_T_start_T_horiz1$nTP + roc_results_auc_T_start_T_horiz1$nFN)
  roc_results_auc_T_start_T_horiz1$fpr = roc_results_auc_T_start_T_horiz1$nFP / (roc_results_auc_T_start_T_horiz1$nTN + roc_results_auc_T_start_T_horiz1$nFP)
  
  #############
  # Then we get prediction error results
  ############
  ape = newdata.id$real_period_status * abs(1 - newdata.id$cum_risk_T_start_T_horiz)  + 
    (1 - newdata.id$real_period_status) * abs(1 - (1 - newdata.id$cum_risk_T_start_T_horiz))
  mape = sum(ape)/nrow(newdata.id)
  
  spe = newdata.id$real_period_status * (1 - newdata.id$cum_risk_T_start_T_horiz)^2  + 
    (1 - newdata.id$real_period_status) * (1 - (1 - newdata.id$cum_risk_T_start_T_horiz))^2
  mspe = sum(spe)/nrow(newdata.id)
  rmspe = sqrt(mspe)
  
  return(list(newdata.id=newdata.id, roc_results=roc_results, 
              roc_results_auc_T_start_T_horiz = roc_results_auc_T_start_T_horiz,
              roc_results_auc_T_start_T_horiz1 = roc_results_auc_T_start_T_horiz1,
              ape = ape, mape=mape, spe = spe,
              mspe = mspe, rmspe = rmspe))
}

environment(goodness_of_fit) = asNamespace("JMbayes")
