#When nb_offset_ratio = 2 it means biopsy distance will be multiplied by 2
minDistScheduleDecision = function(object, patient_data, cur_visit_time=NA,
                                   latest_survival_time=NA, 
                                   nb_offset_ratio = 1, weight_by_horizon_risk=T){

  if(is.na(cur_visit_time)){
    cur_visit_time = min(max(patient_data$year_visit), MAX_FOLLOW_UP)
  }
  
  risk_thresholds = seq(0,1, length.out = 1000)
  res = compareSchedules(object, patient_data, cur_visit_time, latest_survival_time, risk_thresholds, 
                         no_fixed=T, weight_by_horizon_risk = weight_by_horizon_risk)
  
  expected_delays = as.numeric(sapply(res$schedules, "[[", "expected_delay"))
  total_biopsies = as.numeric(sapply(lapply(res$schedules, "[[", "proposed_biopsy_times"), length))
  
  #Distance from optimal point
  dist = sqrt(expected_delays^2 + nb_offset_ratio * (total_biopsies - 1)^2)
  
  #There are many, but we select only one for now
  #optimal_threshold_indices = which(dist == min(dist))
  optimal_threshold_index = which.min(dist)
  next_biopsy_time = res$schedules[[optimal_threshold_index]]$proposed_biopsy_times[1]
  
  return(next_biopsy_time == cur_visit_time)
}


                       
