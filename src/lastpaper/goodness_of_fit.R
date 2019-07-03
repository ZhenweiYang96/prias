goodness_of_fit <- function (object, newdata, T_start, T_horiz, latest_biosy_time_tol = 0.25, horizon = 10, M = 200) {
  
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
  subset_patients_cum_risk_Thoriz = by(INDICES = subset_patients$P_ID, data = subset_patients, FUN = function(pat_data){
    rowMeans(1 - getExpectedFutureOutcomes(object, pat_data, 
                                           latest_survival_time = pat_data$right_cens_time[1],
                                           earliest_failure_time = Inf,
                                           survival_predict_times = T_horiz,
                                           psaDist = "Tdist")$predicted_surv_prob)
  })
  
  newdata.id$real_period_status[newdata.id$reclassification==F & 
                                  newdata.id$right_cens_time < T_horiz] = as.numeric(subset_patients_cum_risk_Thoriz)
  
  #Now for each of the patients obtain their cumulative risk at 
  #various follow-ups under the condition
  #that they didnt have an event until T_start
  cum_risk_T_start_T_horiz = by(INDICES = newdata$P_ID, data = newdata, FUN = function(pat_data){
    rowMeans(1 - getExpectedFutureOutcomes(object, pat_data, 
                                           latest_survival_time = T_start,
                                           earliest_failure_time = Inf,
                                           survival_predict_times = T_horiz,
                                           psaDist = "Tdist")$predicted_surv_prob)
  })
  
  newdata.id$cum_risk_T_start_T_horiz = as.numeric(cum_risk_T_start_T_horiz)
  
  ###############
  # First we get ROC results
  ################
  thresholds = seq(0,1, by=0.001)
  roc_results = t(sapply(thresholds, FUN = function(threshold){
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
              ape = ape, mape=mape, spe = spe,
              mspe = mspe, rmspe = rmspe))
}

environment(goodness_of_fit) = asNamespace("JMbayes")

tt = goodness_of_fit(mvJoint_psa_time_scaled, prias_long_final, 3, 4)