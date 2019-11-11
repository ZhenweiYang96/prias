#When nb1_to_delay = 2 it means 1 biopsy is equal to 2 units of delay
minDistScheduleDecision = function(object, patient_data, cur_visit_time=NA,
                                   latest_survival_time=NA, delay_limit=Inf,
                                   nb1_to_delay = 1, weight_by_horizon_risk=T){
  
  if(is.na(cur_visit_time)){
    cur_visit_time = min(max(patient_data$year_visit), MAX_FOLLOW_UP)
  }
  
  risk_thresholds = seq(from = 0, to = 1, length.out = 1000)
  res = compareSchedules(object, patient_data, cur_visit_time, latest_survival_time, risk_thresholds, 
                         no_fixed=T, weight_by_horizon_risk = weight_by_horizon_risk)
  
  expected_delays = as.numeric(sapply(res$schedules, "[[", "expected_delay"))
  total_biopsies = as.numeric(sapply(lapply(res$schedules, "[[", "proposed_biopsy_times"), length))
  
  #Distance from optimal point
  dist = sqrt(expected_delays^2 + (nb1_to_delay * (total_biopsies - 1))^2)
  
  while(sum(expected_delays<=delay_limit) == 0){
    delay_limit = delay_limit + 0.1
  }
  
  risk_thresholds = risk_thresholds[expected_delays<=delay_limit]
  dist = dist[expected_delays<=delay_limit]
  total_biopsies = total_biopsies[expected_delays<=delay_limit]
  res$schedules = res$schedules[expected_delays<=delay_limit]
  
  #There are many, but we select only one for now
  #optimal_threshold_indices = which(dist == min(dist))
  optimal_threshold_index = which.min(dist)[1]
  next_biopsy_time = res$schedules[[optimal_threshold_index]]$proposed_biopsy_times[1]
  
  return(list(result = (next_biopsy_time == cur_visit_time),
              latest_survival_time=latest_survival_time,
              cur_visit_time = cur_visit_time,
              dist=dist, expected_delays=expected_delays, total_biopsies=total_biopsies,
              risk_thresholds=risk_thresholds))
}