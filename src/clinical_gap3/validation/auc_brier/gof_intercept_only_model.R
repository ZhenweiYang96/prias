goodness_of_fit_intercept_only <- function (icenRegModel, newdata, T_start, T_horiz) {
  
  #we want intercept only model estimates
  baseline_est = getSCurves(icenRegModel, newdata=NULL)
  baseline_est_int = baseline_est$Tbull_ints
  baseline_est_int = cbind(baseline_est_int, baseline_est$S_curves$baseline)
  baseline_est_int = rbind(c(0,0,1), baseline_est_int)
  
  getConditionalSurvProb = function(last.time, pred.time){
    if(pred.time<=last.time){
      return(1)
    }else{
      #turnbull intervals are hard to understand
      index1 = tail(which(baseline_est_int[,2]<=last.time),1)
      if(baseline_est_int[index1 + 1, 1]<=last.time){
        last.time.surv = baseline_est_int[index1 + 1,3]
      }else{
        last.time.surv = baseline_est_int[index1, 3]        
      }
      
      index2 = tail(which(baseline_est_int[,2]<=pred.time),1)
      if(baseline_est_int[index2 + 1, 1]<=pred.time){
        pred.time.surv = baseline_est_int[index2 + 1,3]
      }else{
        pred.time.surv = baseline_est_int[index2, 3]
      }
      
      return(pred.time.surv / last.time.surv)
    }
  }
  
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
        newdata.id$real_period_status[newdata.id$P_ID == p_id] = 1 - getConditionalSurvProb(last.time = pat_data$right_cens_time[1], pred.time = T_horiz)
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
      newdata.id$cum_risk_T_start_T_horiz[newdata.id$P_ID == p_id] = 1 - getConditionalSurvProb(last.time = T_start, pred.time = T_horiz)
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