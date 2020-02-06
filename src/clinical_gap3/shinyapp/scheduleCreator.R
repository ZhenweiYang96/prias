getFixedSchedule = function(latest_survival_time=NA,
                            min_biopsy_gap = 1, biopsy_frequency=1){
  
  visit_schedule = patient_cache$visit_schedule
  
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

getPRIASSchedule = function(latest_survival_time, min_biopsy_gap){
  
  visit_schedule = patient_cache$visit_schedule
  psa_indices = sapply(visit_schedule, FUN = function(x){
    which(patient_cache$PSA_CACHE_TIMES==x)
  })
  
  future_log2psaplus1_matrix = patient_cache$PSA_CACHE_FULL[psa_indices,,drop = FALSE]
    
  proposed_biopsy_times_list = vector("list", length=max(min(100, M)))
  for(sim in 1:length(proposed_biopsy_times_list)){
    proposed_biopsy_times = c()
    
    future_log2psaplus1 = apply(future_log2psaplus1_matrix, MARGIN = 1, sample, size=1)
    future_psa = pmax(0.1, 2^(future_log2psaplus1)-1)
    obs_psa = patient_cache$patient_data$psa
    obs_psa_times = patient_cache$patient_data$year_visit
    
    if(max(obs_psa_times)==patient_cache$current_visit_time){
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

getRiskBasedSchedule = function(surv_threshold, previous_test_time, gap){
  proposed_test_times <- c()
  surv_schedule_temp <- patient_cache$SURV_CACHE_FULL
  fixed_grid_visits <- patient_cache$visit_schedule
  
  for(i in 1:length(fixed_grid_visits)){
    if(fixed_grid_visits[i] - previous_test_time>=gap){
      surv_cache_index = which(patient_cache$SURV_CACHE_TIMES==fixed_grid_visits[i])
      surv_row_visit = surv_schedule_temp[surv_cache_index,]
      if(mean(surv_row_visit, na.rm = T) <= surv_threshold){
        proposed_test_times <- c(proposed_test_times, fixed_grid_visits[i])
        previous_test_time = fixed_grid_visits[i]
        surv_schedule_temp <- apply(patient_cache$SURV_CACHE_FULL, 2, FUN = function(x){x/x[surv_cache_index]})
      }
    }
  }
  
  if(length(proposed_test_times)==0){
    proposed_test_times = MAX_FOLLOW_UP
  }
  
  return(proposed_test_times)
}

getAllRiskSchedule = function(previous_test_time, gap){
  surv_thresholds = seq(0,1,by=0.01)
  
  proposed_test_times = lapply(surv_thresholds, FUN = getRiskBasedSchedule, previous_test_time, gap)
  
  proposed_test_time_str = sapply(proposed_test_times, paste0, collapse="-")
  unique_test_time_str = unique(proposed_test_time_str)
  unique_test_time_list = lapply(strsplit(unique(sapply(proposed_test_times, paste0, collapse="-")), split = "-"), as.numeric)
  unique_consequences = lapply(unique_test_time_list, getConsequences, gap, previous_test_time)
  names(unique_consequences) = unique_test_time_str
  
  total_consequences = lapply(proposed_test_time_str, FUN = function(str){
    unique_consequences[[str]]
  })
  names(total_consequences) = surv_thresholds
  
  return(total_consequences)
}  

getConsequences = function(proposed_test_times, gap, last_test_time){
  set.seed(2019);
  
  if(length(proposed_test_times)==0){
    proposed_test_times = MAX_FOLLOW_UP
  }
  
  planned_test_schedule = proposed_test_times
  #Make the last test at MAX_FOLLOW_UP
  if(MAX_FOLLOW_UP - tail(planned_test_schedule,1) < gap){
    planned_test_schedule = planned_test_schedule[-length(planned_test_schedule)]
  }
  planned_test_schedule = c(planned_test_schedule, MAX_FOLLOW_UP)
  
  #Now we create intervals in which to integrate the conditional surv prob
  if(length(planned_test_schedule)==1){
    test_intervals = list(c(last_test_time, planned_test_schedule))
  }else{
    test_intervals = c(last_test_time,
                       rep(planned_test_schedule[-length(planned_test_schedule)],each=2),
                       planned_test_schedule[length(planned_test_schedule)])
    test_intervals = split(test_intervals, rep(1:(length(test_intervals)/2), each=2))
  }
  
  SURV_CACHE_FULL_rescaled = 1 - t(apply(1-patient_cache$SURV_CACHE_FULL, 1, FUN = function(x){
    x / (1-patient_cache$SURV_CACHE_FULL[nrow(patient_cache$SURV_CACHE_FULL),, drop=F])
  })) 
  
  interval_res = vector("list", length(test_intervals))
  exp_num_tests = rep(0, M)
  for(j in 1:length(test_intervals)){
    wt_points=getGaussianQuadWeightsPoints(test_intervals[[j]])
    lower_limit = test_intervals[[j]][1]
    upper_limit = test_intervals[[j]][2]
    
    lower_limit_nearest_index = which.min(abs(lower_limit - patient_cache$SURV_CACHE_TIMES))
    upper_limit_nearest_index = which.min(abs(upper_limit - patient_cache$SURV_CACHE_TIMES))
    interval_res[[j]]$cum_risk_interval = SURV_CACHE_FULL_rescaled[lower_limit_nearest_index,] - 
      SURV_CACHE_FULL_rescaled[upper_limit_nearest_index,]
    
    cond_expected_fail_time = sapply(1:length(wt_points$points), function(i){
      cum_surv_at_points = SURV_CACHE_FULL_rescaled[which.min(abs(wt_points$points[i]-patient_cache$SURV_CACHE_TIMES)),] - SURV_CACHE_FULL_rescaled[upper_limit_nearest_index,]
      scaled_cum_surv_at_points = cum_surv_at_points/interval_res[[j]]$cum_risk_interval
      
      return(wt_points$weights[i] * scaled_cum_surv_at_points)
    })
    
    interval_res[[j]]$delay = upper_limit - (lower_limit + apply(cond_expected_fail_time, 1, sum))
    
    exp_num_tests = exp_num_tests + j * interval_res[[j]]$cum_risk_interval
  }
  exp_num_tests = exp_num_tests + j * SURV_CACHE_FULL_rescaled[upper_limit_nearest_index,]
  exp_num_tests = mean(exp_num_tests, na.rm = T)
  
  expected_detection_delay = mean(apply(sapply(interval_res, FUN = function(x){
    x$delay * x$cum_risk_interval
  }),1, sum, na.rm=T),na.rm=T)
  
  #proposed test times always has a final test at MAX_FOLLOW_UP
  return(list(expected_detection_delay = expected_detection_delay,
              expected_num_tests = exp_num_tests,
              planned_test_schedule = planned_test_schedule))
}
