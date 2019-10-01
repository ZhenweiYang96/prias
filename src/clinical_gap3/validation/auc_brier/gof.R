goodness_of_fit <- function (pred_model, orig_model, newdata, T_start, T_horiz, M = 500) {
  
  orig_data = newdata
  
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
  newdata.id.subset = droplevels(newdata.id[newdata.id$reclassification==F & newdata.id$right_cens_time < T_horiz,])
  
  if(nrow(newdata.id.subset) > 0){
    for(p_id in newdata.id.subset$P_ID){
      pat_data = orig_data[orig_data$P_ID == p_id,]
      if(all(is.na(pat_data$log2psaplus1))){
        newdata.id$real_period_status[newdata.id$P_ID == p_id] = NA
      }else{
        newdata.id$real_period_status[newdata.id$P_ID == p_id] = rowMeans(1 - getExpectedFutureOutcomes(orig_model, pat_data, 
                                                                                                        latest_survival_time = pat_data$right_cens_time[1],
                                                                                                        earliest_failure_time = Inf,
                                                                                                        survival_predict_times = T_horiz,
                                                                                                        psaDist = "Tdist", M=M)$predicted_surv_prob, na.rm = T)
      }
    }
  }
  
  #various follow-ups under the condition
  #that they didnt have an event until T_start
  for(p_id in newdata.id$P_ID){
    pat_data = newdata[newdata$P_ID == p_id,]
    
    if(all(is.na(pat_data$log2psaplus1))){
      newdata.id$cum_risk_T_start_T_horiz[newdata.id$P_ID == p_id] = NA
    }else{
      newdata.id$cum_risk_T_start_T_horiz[newdata.id$P_ID == p_id] = rowMeans(1 - getExpectedFutureOutcomes(pred_model, pat_data, 
                                                                                                            latest_survival_time = T_start,
                                                                                                            earliest_failure_time = Inf,
                                                                                                            survival_predict_times = T_horiz,
                                                                                                            psaDist = "Tdist", M = M)$predicted_surv_prob, na.rm = T)
    }
  }
  
  #############
  # Then we get prediction error results
  ############
  newdata.id$ape = newdata.id$real_period_status * abs(1 - newdata.id$cum_risk_T_start_T_horiz)  + 
    (1 - newdata.id$real_period_status) * abs(1 - (1 - newdata.id$cum_risk_T_start_T_horiz))
  
  newdata.id$spe = newdata.id$real_period_status * (1 - newdata.id$cum_risk_T_start_T_horiz)^2  + 
    (1 - newdata.id$real_period_status) * (1 - (1 - newdata.id$cum_risk_T_start_T_horiz))^2
  
  return(newdata.id=newdata.id)
}

environment(goodness_of_fit) = asNamespace("JMbayes")
