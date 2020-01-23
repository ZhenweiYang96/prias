#These schedules do not check if minimum gap is maintained between biopsies
ifFixedBiopsy = function(cur_visit_time, last_biopsy_time,
                         biopsy_frequency = 1){
  return(cur_visit_time - last_biopsy_time>=biopsy_frequency)
}

getFixedSchedule = function(cur_visit_time, biopsy_frequency=1, horizon=10){
  return(seq(cur_visit_time, horizon, by=biopsy_frequency))
}

#here we check if minimum gap is maintained
ifPRIASBiopsy = function(patient_data, cur_visit_time, last_biopsy_time){
  
  if(cur_visit_time - last_biopsy_time < 1){
    return(FALSE)
  }else{
    fixed_schedule = c(1, 4, 7, 10, 15)
    
    if(cur_visit_time %in% fixed_schedule){
      return(TRUE)
    }else{
      psa = pmax(0.1, 2^(patient_data$log2psaplus1)-1)
      log2psa = log(psa, base = 2)
      year_visit = patient_data$year_visit
      
      psa_dt = 1/(lm(log2psa~year_visit)$coefficients[2])
      
      return(psa_dt>=0 & psa_dt<=10)
    }
  }
}  

#assuming that observed longitudinal data is available up to current visit time
getPRIASSchedule = function(object, patient_data, cur_visit_time, 
                            last_biopsy_time,
                            M=500, horizon=10){
  
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, 10, 0.5))
  visit_schedule = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME <=horizon])
  future_log2psaplus1_matrix = getExpectedFutureOutcomes(object, patient_data, last_biopsy_time, 
                                                         long_predict_times = visit_schedule, 
                                                         psaDist = "Tdist", M = M, addRandomError = T)$predicted_psa
  future_log2psaplus1_matrix[1,] = rep(patient_data$log2psaplus1[nrow(patient_data)], M)
  patient_data = patient_data[-nrow(patient_data)]
  
  proposed_biopsy_times_list = vector("list", length=500)
  for(sim in 1:length(proposed_biopsy_times_list)){
    proposed_biopsy_times = c()
    temp_patient_data = patient_data
    previous_biopsy_time = last_biopsy_time
    
    for(i in 1:length(visit_schedule)){
      new_row = temp_patient_data[1,]
      new_row$year_visit = visit_schedule[i]
      new_row$log2psaplus1 = sample(x=future_log2psaplus1_matrix[i,], size=1)
      
      temp_patient_data = rbind(temp_patient_data, new_row)
      
      if(ifPRIASBiopsy(temp_patient_data, visit_schedule[i], previous_biopsy_time)){
        proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
        previous_biopsy_time = visit_schedule[i]
      }
    }
    
    proposed_biopsy_times_list[[sim]] = proposed_biopsy_times
  }
  
  mode_schedule_string = names(table(unlist(lapply(proposed_biopsy_times_list, paste, collapse = "-"))))[1]
  
  final_proposed_biopsy_times = as.numeric(strsplit(mode_schedule_string, split = '-')[[1]])
  return(final_proposed_biopsy_times)
}

ifFixedRiskBasedBiopsy =  function(object, patient_data, cur_visit_time, 
                                   last_biopsy_time, min_biopsy_gap=1, threshold, M=500){
  
  if(cur_visit_time - last_biopsy_time>=min_biopsy_gap){
    visit_cum_risk = rowMeans(1-getExpectedFutureOutcomes(object, patient_data, last_biopsy_time, 
                                                          survival_predict_times = cur_visit_time,
                                                          psaDist = "Tdist", M = M)$predicted_surv_prob)
    return(visit_cum_risk >= threshold)
  }else{
    return(FALSE)
  }
}

getFixedRiskBasedScheduleNoGrid = function(object, patient_data, 
                                           cur_visit_time, 
                                           last_biopsy_time,
                                           risk_threshold = 0.05, 
                                           min_biopsy_gap = 1,
                                           M=500, horizon=10){
  previous_biopsy_time = last_biopsy_time
  visit_time = cur_visit_time
  proposed_biopsy_times = c()
  while(horizon - previous_biopsy_time >= min_biopsy_gap){
    if(visit_time - previous_biopsy_time >= min_biopsy_gap){
      grid_times = seq(visit_time, horizon, by=0.01)
      visit_cum_risk = rowMeans(1-getExpectedFutureOutcomes(object, patient_data, previous_biopsy_time, 
                                                            survival_predict_times = grid_times,
                                                            psaDist = "Tdist", M = M)$predicted_surv_prob)
      
      nearest_index = which.min(abs(visit_cum_risk - risk_threshold))
      previous_biopsy_time = grid_times[nearest_index]
      
      proposed_biopsy_times = c(proposed_biopsy_times, previous_biopsy_time)
    }
    
    visit_time = previous_biopsy_time + min_biopsy_gap
  }
  
  return(proposed_biopsy_times)
}

#################### CACHE BASED SOLUTION #######################
ifAutomaticRiskBasedBiopsy = function(object, patient_data, cur_visit_time,
                                      last_biopsy_time, min_biopsy_gap = 1, M = M, 
                                      delay_limit = 1, horizon=10){
  if(cur_visit_time - last_biopsy_time < min_biopsy_gap){
    return(FALSE)
  }
  
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
  
  #Step 1: Create a CACHE of predicted survival probabilities to speed up operation
  if(horizon - last_biopsy_time <= 1.5){
    SURV_CACHE_TIMES = c(seq(last_biopsy_time, horizon, by = 1/365), horizon)
  }else{
    CACHE_SIZE = 500
    SURV_CACHE_TIMES = seq(last_biopsy_time, horizon, length.out = CACHE_SIZE)
  }
  SURV_CACHE_TIMES <- unique(sort(c(PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME >= last_biopsy_time], SURV_CACHE_TIMES), decreasing = F))
  
  pred_res = getExpectedFutureOutcomes(object, patient_data, last_biopsy_time, 
                                       survival_predict_times = SURV_CACHE_TIMES[-1],
                                       psaDist = "Tdist", M = M)
  SURV_CACHE_FULL = rbind(rep(1, M), pred_res$predicted_surv_prob)
  rownames(SURV_CACHE_FULL)[1] = as.character(last_biopsy_time)
  
  fixed_grid_visits = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME<=horizon])
  
  #Fixed risk schedule function
  getRiskScheduleWithCache <- function(surv_threshold){
    proposed_biopsy_times <- c()
    previous_biopsy_time <- last_biopsy_time
    surv_schedule_temp <- SURV_CACHE_FULL
    
    for(i in 1:length(fixed_grid_visits)){
      if(fixed_grid_visits[i] - previous_biopsy_time>=min_biopsy_gap){
        surv_cache_index = which(SURV_CACHE_TIMES==fixed_grid_visits[i])
        surv_row_visit = surv_schedule_temp[surv_cache_index,]
        if(mean(surv_row_visit, na.rm = T) <= surv_threshold){
          proposed_biopsy_times <- c(proposed_biopsy_times, fixed_grid_visits[i])
          previous_biopsy_time = fixed_grid_visits[i]
          surv_schedule_temp <- apply(SURV_CACHE_FULL, 2, FUN = function(x){x/x[surv_cache_index]})
        }
      }
    }
    
    return(proposed_biopsy_times)
  }
  
  risk_thresholds = seq(0, 1, length.out = 201)
  res = vector("list", length=length(risk_thresholds))
  names(res) = c(paste("Risk:", risk_thresholds))
  
  for(i in 1:length(risk_thresholds)){
    threshold = risk_thresholds[i]
    risk_schedule = getRiskScheduleWithCache(1-threshold)
    if(is.null(risk_schedule)){
      risk_schedule = horizon
    }
    
    res[[i]] = getConsequencesWithCache(SURV_CACHE_TIMES, SURV_CACHE_FULL, horizon,
                                        risk_schedule, last_biopsy_time)
    res[[i]]$threshold = threshold
  }
  
  expected_delays = sapply(res, "[[", "expected_delay")
  exp_total_biopsies = sapply(res, "[[", "expected_num_biopsies")
  
  #Distance from optimal point
  dist = sqrt(expected_delays^2 + (exp_total_biopsies - 1)^2)
  
  while(sum(expected_delays<=delay_limit) == 0){
    delay_limit = delay_limit + 0.1
  }
  
  filter = expected_delays<=delay_limit
  risk_thresholds = risk_thresholds[filter]
  dist = dist[filter]
  exp_total_biopsies = exp_total_biopsies[filter]
  res = res[filter]
  expected_delays = expected_delays[filter]
  
  #There are many, but we select only one for now
  optimal_schedule_index = which.min(dist)[1]
  optimal_schedule = res[[optimal_schedule_index]]$practical_biopsy_times
  decision = optimal_schedule[1] <= cur_visit_time
  
  return(decision)
}

ifAutomaticFullTreeBasedBiopsy = function(object, patient_data, cur_visit_time,
                                          last_biopsy_time, min_biopsy_gap = 1, M = M, 
                                          delay_limit = 1, horizon=10){
  if(cur_visit_time - last_biopsy_time < min_biopsy_gap){
    return(FALSE)
  }
  
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
  
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
  
  #Step 1: Create a CACHE of predicted survival probabilities to speed up operation
  if(horizon - last_biopsy_time <= 1.5){
    SURV_CACHE_TIMES = c(seq(last_biopsy_time, horizon, by = 1/365), horizon)
  }else{
    CACHE_SIZE = 500
    SURV_CACHE_TIMES = seq(last_biopsy_time, horizon, length.out = CACHE_SIZE)
  }
  SURV_CACHE_TIMES <- unique(sort(c(PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME >= last_biopsy_time], SURV_CACHE_TIMES), decreasing = F))
  
  pred_res = getExpectedFutureOutcomes(object, patient_data, last_biopsy_time, 
                                       survival_predict_times = SURV_CACHE_TIMES[-1],
                                       psaDist = "Tdist", M = M)
  
  SURV_CACHE_FULL = rbind(rep(1, M), pred_res$predicted_surv_prob)
  
  #Function to create all possible biopsy schedules 
  allBiopsySchedules = function(cur_visit_time, min_biopsy_gap, horizon){
    times = list()
    count = 0
    
    PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
    
    recursive = function(last_biopsy_time, current_visit_time, cur_path){
      if(current_visit_time < cur_visit_time){
        count <<- count + 1
        times[[count]] <<- as.numeric(strsplit(cur_path, split = "-")[[1]])
      }else{
        if(last_biopsy_time - current_visit_time>=min_biopsy_gap){
          #case 1 do a biopsy
          recursive(current_visit_time, 
                    current_visit_time - min_biopsy_gap,
                    paste0(current_visit_time,"-",cur_path))
          
          #case 2 don't do a biopsy
          recursive(last_biopsy_time, 
                    tail(PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME<current_visit_time],1),
                    cur_path)
        }
      }
    }
    recursive(horizon, horizon-min_biopsy_gap, as.character(horizon))
    
    #I will sort them by first biopsy time
    first_biopsy_time = sapply(times, min)
    total_biopsies = sapply(times, length)
    times = times[order(first_biopsy_time, total_biopsies, decreasing = F)]
    return(times)
  }
  
  allSchedules = allBiopsySchedules(cur_visit_time = cur_visit_time, 
                                    min_biopsy_gap = min_biopsy_gap, horizon = horizon)
  
  res = vector("list", length=length(allSchedules))
  for(i in 1:length(allSchedules)){
    schedule = allSchedules[[i]]
    
    res[[i]] = getConsequencesWithCache(SURV_CACHE_TIMES, SURV_CACHE_FULL, horizon,
                                        schedule, last_biopsy_time)
  }
  
  expected_delays = sapply(res, "[[", "expected_delay")
  exp_total_biopsies = sapply(res, "[[", "expected_num_biopsies")
  
  #Distance from optimal point
  dist = sqrt(expected_delays^2 + (exp_total_biopsies - 1)^2)
  
  while(sum(expected_delays<=delay_limit) == 0){
    delay_limit = delay_limit + 0.1
  }
  
  risk_thresholds = risk_thresholds[expected_delays<=delay_limit]
  dist = dist[expected_delays<=delay_limit]
  exp_total_biopsies = exp_total_biopsies[expected_delays<=delay_limit]
  res = res[expected_delays<=delay_limit]
  
  #There are many, but we select only one for now
  optimal_schedule_index = which.min(dist)[1]
  optimal_schedule = res[[optimal_schedule_index]]$practical_biopsy_times
  decision = optimal_schedule[1] <= cur_visit_time
  
  return(decision)
}

getConsequencesWithCache = function(cache_times, cache_surv_full, horizon,
                                    proposed_biopsy_times, latest_survival_time){
  
  practical_biopsy_times = proposed_biopsy_times
  #Make the last biopsy at year 10
  if(horizon - tail(practical_biopsy_times,1) < 1){
    practical_biopsy_times = practical_biopsy_times[-length(practical_biopsy_times)]
  }
  practical_biopsy_times = c(practical_biopsy_times, horizon)
  
  #Now we create intervals in which to integrate the conditional surv prob
  if(length(practical_biopsy_times)==1){
    biop_intervals = list(c(latest_survival_time, practical_biopsy_times))
  }else{
    biop_intervals = c(latest_survival_time,
                       rep(practical_biopsy_times[-length(practical_biopsy_times)],each=2),
                       practical_biopsy_times[length(practical_biopsy_times)])
    biop_intervals = split(biop_intervals, rep(1:(length(biop_intervals)/2), each=2))
  }
  
  res = vector("list", length(biop_intervals))
  exp_num_biopsies = rep(0, M)
  for(j in 1:length(biop_intervals)){
    wt_points=getGaussianQuadWeightsPoints(biop_intervals[[j]])
    lower_limit = biop_intervals[[j]][1]
    upper_limit = biop_intervals[[j]][2]
    
    lower_limit_nearest_index = which(lower_limit==cache_times)
    upper_limit_nearest_index = which(upper_limit==cache_times)
    res[[j]]$cum_risk_interval = cache_surv_full[lower_limit_nearest_index,] - 
      cache_surv_full[upper_limit_nearest_index,]
    
    cond_expected_fail_time = sapply(1:length(wt_points$points), function(i){
      cum_surv_at_points = cache_surv_full[which.min(abs(wt_points$points[i]-cache_times)),] - cache_surv_full[upper_limit_nearest_index,]
      scaled_cum_surv_at_points = cum_surv_at_points/res[[j]]$cum_risk_interval
      
      return(wt_points$weights[i] * scaled_cum_surv_at_points)
    })
    
    res[[j]]$delay = upper_limit - (lower_limit + apply(cond_expected_fail_time, 1, sum))
    exp_num_biopsies = exp_num_biopsies + j * res[[j]]$cum_risk_interval
  }
  exp_num_biopsies = exp_num_biopsies + j * cache_surv_full[upper_limit_nearest_index,]
  exp_num_biopsies = mean(exp_num_biopsies, na.rm = T)
  
  expected_delay = mean(apply(sapply(res, FUN = function(x){
    x$delay * x$cum_risk_interval
  }),1, sum, na.rm=T),na.rm=T)
  
  #proposed biopsy times always has a final biopsy at year 10
  return(list(expected_delay = expected_delay,
              expected_num_biopsies = exp_num_biopsies,
              proposed_biopsy_times=proposed_biopsy_times,
              practical_biopsy_times = practical_biopsy_times))
}