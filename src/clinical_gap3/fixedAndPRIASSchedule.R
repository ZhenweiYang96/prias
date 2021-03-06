getFixedSchedule = function(cur_visit_time, latest_survival_time=NA,
                            min_biopsy_gap = 1, biopsy_frequency=1, horizon=10){
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, 10, 0.5))
  
  visit_schedule = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME <=horizon])
  
  proposed_biopsy_times = c()
  min_biopsy_gap = max(min_biopsy_gap, biopsy_frequency)
  
  latest_biopsy_time = latest_survival_time
  for(i in 1:length(visit_schedule)){
    visit_time = visit_schedule[i]
    if(visit_time - latest_biopsy_time>=min_biopsy_gap){
      proposed_biopsy_times = c(proposed_biopsy_times, visit_time)
      latest_biopsy_time = visit_time
    }
  }
  
  return(proposed_biopsy_times)
}

getPRIASSchedule = function(object, patient_data, cur_visit_time=NA, 
                            latest_survival_time=NA,
                            min_biopsy_gap = 1,
                            M=750, horizon=10){
  if(is.na(cur_visit_time)){
    cur_visit_time = min(max(patient_data$year_visit), horizon)
  }
  
  patient_data = patient_data[patient_data$year_visit <= cur_visit_time,]
  
  if(is.na(latest_survival_time)){
    latest_survival_time = max(patient_data$year_visit[!is.na(patient_data$gleason_sum)])
  }
  
  if(cur_visit_time < latest_survival_time){
    stop("Current visit time should be more than latest survival time")
  }
  
  #making schedule now
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, 10, 0.5))
  visit_schedule = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME <=horizon])
  future_log2psaplus1_matrix = getExpectedFutureOutcomes(object, patient_data, latest_survival_time, 
                                                         psa_predict_times = visit_schedule, 
                                                         psaDist = "Tdist", M = M, addRandomError = T)$predicted_psa
  proposed_biopsy_times_list = vector("list", length=500)
  for(sim in 1:length(proposed_biopsy_times_list)){
    proposed_biopsy_times = c()
    
    future_log2psaplus1 = apply(future_log2psaplus1_matrix, MARGIN = 1, sample, size=1)
    future_psa = pmax(0.1, 2^(future_log2psaplus1)-1)
    obs_psa = patient_data$psa
    obs_psa_times = patient_data$year_visit
    
    if(max(obs_psa_times)==cur_visit_time){
      future_psa[1] = tail(obs_psa,1)
      obs_psa = obs_psa[-length(obs_psa)]
      obs_psa_times = obs_psa_times[-length(obs_psa_times)]
    }
    
    fixed_schedule = c(1, 4, 7, 10, 15)
    latest_biopsy_time = latest_survival_time
    for(i in 1:length(visit_schedule)){
      
      if(latest_biopsy_time >= fixed_schedule[1]){
        fixed_schedule = fixed_schedule[-1]
      }
      
      #in some cases two PSA at current visit time 
      log2psa = log(c(obs_psa, future_psa[1:i]), base = 2)
      year_visit = c(obs_psa_times, visit_schedule[1:i])
      
      psa_dt = 1/(lm(log2psa~year_visit)$coefficients[2])
      
      if((visit_schedule[i] - latest_biopsy_time) >= min_biopsy_gap){
        #if switch to annual schedule
        if(psa_dt>=0 & psa_dt<=10){
          latest_biopsy_time = visit_schedule[i]
          proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
        } else if(visit_schedule[i] >= fixed_schedule[1]){
          latest_biopsy_time = visit_schedule[i]
          proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
        }
      }
    }
    proposed_biopsy_times_list[[sim]] = proposed_biopsy_times
  }
  
  mode_schedule_string = names(table(unlist(lapply(proposed_biopsy_times_list, paste, collapse = "-"))))[1]
  
  final_proposed_biopsy_times = as.numeric(strsplit(mode_schedule_string, split = '-')[[1]])
  return(final_proposed_biopsy_times)
}