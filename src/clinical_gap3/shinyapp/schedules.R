getPRIASSchedule = function(latest_survival_time, obs_psa, obs_psa_times, 
                            visit_schedule, future_log2psaplus1){
  
  fixed_schedule = c(1,4,7,10,15)
  
  proposed_biopsy_times = c()
  latest_biopsy_time = latest_survival_time
  
  for(i in 1:length(visit_schedule)){
    log2psa = log(c(obs_psa, 2^(future_log2psaplus1[1:i])-1), base = 2)
    year_visit = c(obs_psa_times,visit_schedule[1:i])
    
    psa_dt = 1/(lm(log2psa~year_visit)$coefficients[2])
    
    #if switch to annual schedule
    if(psa_dt>=0 & psa_dt<=10){
      if((visit_schedule[i] - latest_biopsy_time) >= 1){
        latest_biopsy_time = visit_schedule[i]
        proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
      }
      #else wait
    }else{
      if(visit_schedule[i] %in% fixed_schedule){
        if((visit_schedule[i] - latest_biopsy_time) >= 1){
          latest_biopsy_time = visit_schedule[i]
          proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
        }
        #else wait
      }
    }
  }
  
  return(proposed_biopsy_times)
}

getRiskSchedule = function(surv_threshold,
                           latest_survival_time,
                           visit_schedule, surv_schedule,
                           min_biopsy_gap){
  proposed_biopsy_times = c()
  latest_biopsy_time = latest_survival_time
  
  surv_schedule_temp = surv_schedule
  for(i in 1:length(visit_schedule)){
    if(mean(surv_schedule_temp[i,]) <= surv_threshold & 
       (visit_schedule[i]-latest_biopsy_time)>=min_biopsy_gap){
      latest_biopsy_time = visit_schedule[i]
      proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
      surv_schedule_temp = t(apply(surv_schedule, 1, FUN = function(x){
        x/surv_schedule[i,]
      }))
    }
  }
  
  return(proposed_biopsy_times)
}

getConsequences = function(cache_times, cache_surv_full, horizon,
                           proposed_biopsy_times, latest_survival_time){
  
  orig_biopsy_times = proposed_biopsy_times
  #Make the last biopsy at year 10
  #and if the last biopsy was 9.5 or more, remove it and only keep the one at 10
  if(max(proposed_biopsy_times)>=(horizon-0.5)){
    proposed_biopsy_times = proposed_biopsy_times[-length(proposed_biopsy_times)]
  }
  proposed_biopsy_times = c(proposed_biopsy_times, horizon)
  
  
  #Now we create intervals in which to integrate the conditional surv prob
  if(length(proposed_biopsy_times)==1){
    biop_intervals = list(c(latest_survival_time,proposed_biopsy_times))
  }else{
    biop_intervals = c(latest_survival_time,
                       rep(proposed_biopsy_times[-length(proposed_biopsy_times)],each=2),
                       proposed_biopsy_times[length(proposed_biopsy_times)])
    biop_intervals = split(biop_intervals, rep(1:(length(biop_intervals)/2), each=2))
  }
  
  res = vector("list", length(biop_intervals))
  for(j in 1:length(biop_intervals)){
    wt_points=getGaussianQuadWeightsPoints(biop_intervals[[j]])
    lower_limit = biop_intervals[[j]][1]
    upper_limit = biop_intervals[[j]][2]
    
    lower_limit_nearest_index = which.min(abs(lower_limit-cache_times))
    upper_limit_nearest_index = which.min(abs(upper_limit-cache_times))
    res[[j]]$cum_risk_interval = cache_surv_full[lower_limit_nearest_index,] - 
      cache_surv_full[upper_limit_nearest_index,]
    
    cond_expected_fail_time = sapply(1:length(wt_points$points), function(i){
      cum_surv_at_points = cache_surv_full[which.min(abs(wt_points$points[i]-cache_times)),] - cache_surv_full[upper_limit_nearest_index,]
      scaled_cum_surv_at_points = cum_surv_at_points/res[[j]]$cum_risk_interval
      
      return(wt_points$weights[i] * scaled_cum_surv_at_points)
    })
    
    res[[j]]$delay = upper_limit - (lower_limit + apply(cond_expected_fail_time,1,sum))
  }
  
  expected_delay = mean(apply(sapply(res, FUN = function(x){
    x$delay * x$cum_risk_interval
  }),1, sum, na.rm=T),na.rm=T)
  
  return(list(expected_delay = expected_delay,
              total_biopsies = length(orig_biopsy_times),
              biopsy_times=orig_biopsy_times))
}
