compareSchedules = function(object, patient_data, cur_visit_time=NA, latest_survival_time=NA,
                            risk_thresholds = c(0.05, 0.1, 0.15), 
                            weight_by_horizon_risk=T,
                            min_biopsy_gap = 1,
                            M=500, CACHE_SIZE=200, no_fixed=F){
  
  if(is.na(cur_visit_time)){
    cur_visit_time = min(max(patient_data$year_visit), MAX_FOLLOW_UP)
  }
  
  patient_data = patient_data[patient_data$year_visit <= cur_visit_time,]
  
  if(is.na(latest_survival_time)){
    latest_survival_time = max(patient_data$year_visit[!is.na(patient_data$gleason_sum)])
  }
  
  if(cur_visit_time < latest_survival_time){
    stop("Current visit time should be more than latest survival time")
  }
  
  
  SURV_CACHE_TIMES = c(latest_survival_time, seq(latest_survival_time+1/365, MAX_FOLLOW_UP, length.out = CACHE_SIZE-1))
  pred_res=getExpectedFutureOutcomes(object, patient_data, latest_survival_time, 
                                     survival_predict_times = SURV_CACHE_TIMES[-1],
                                     psa_predict_times = SURV_CACHE_TIMES, 
                                     psaDist = "Tdist", M = M, addRandomError = T)
  
  PSA_CACHE_FULL = pred_res$predicted_psa
  SURV_CACHE_FULL = rbind(rep(1, M), pred_res$predicted_surv_prob)
  rm(pred_res)
  
  if(weight_by_horizon_risk==F){
    SURV_CACHE_FULL = 1 - t(apply(1-SURV_CACHE_FULL, 1, FUN = function(x){
      x / (1-SURV_CACHE_FULL[CACHE_SIZE,])
    }))
  }
  
  createSchedules(patient_data, cur_visit_time,
                  latest_survival_time, risk_thresholds,
                  min_biopsy_gap, SURV_CACHE_TIMES, SURV_CACHE_TIMES,
                  SURV_CACHE_FULL, PSA_CACHE_FULL, no_fixed)
}

createSchedules = function(patient_data, cur_visit_time, latest_survival_time,
                           risk_thresholds = c(0.05, 0.1, 0.15), 
                           min_biopsy_gap, SURV_CACHE_TIMES, PSA_CACHE_TIMES,
                           SURV_CACHE_FULL, PSA_CACHE_FULL, no_fixed=F){
  
  #Decision epochs in years
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP, 0.5))
  
  horizon_surv_prob = min(rowMeans(SURV_CACHE_FULL, na.rm = T))
  
  visit_schedule = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME<=MAX_FOLLOW_UP])
  visit_nearest_indices = sapply(visit_schedule, FUN = function(x){
    which.min(abs(x-SURV_CACHE_TIMES))
  })
  surv_schedule = SURV_CACHE_FULL[visit_nearest_indices,, drop=F]
  PSA_CACHE_FULL = PSA_CACHE_FULL[visit_nearest_indices,, drop=F]
  
  #added 3 for prias, annual and biennial schedule
  res = vector("list", length=length(risk_thresholds))
  names(res) = c(paste("Risk:", risk_thresholds))
  
  for(i in 1:length(risk_thresholds)){
    threshold = risk_thresholds[i]
    risk_schedule = getRiskSchedule(1-threshold, latest_survival_time, visit_schedule, surv_schedule, min_biopsy_gap)
    res[[i]] = getConsequences(SURV_CACHE_TIMES, SURV_CACHE_FULL, MAX_FOLLOW_UP,
                               risk_schedule, latest_survival_time)
    res[[i]]$threshold = threshold
  }
  
  if(no_fixed==F){
    #first make the annual schedule
    annual_schedule = if(cur_visit_time - latest_survival_time < min_biopsy_gap){
      seq(latest_survival_time + min_biopsy_gap, MAX_FOLLOW_UP, by=1)
    }else{
      seq(cur_visit_time, MAX_FOLLOW_UP, by=1)
    }
    res$annual = getConsequences(SURV_CACHE_TIMES, SURV_CACHE_FULL,MAX_FOLLOW_UP,
                                 annual_schedule, latest_survival_time)
    
    #now make the biennial schedule
    biennial_schedule = if(cur_visit_time - latest_survival_time < 2){
      seq(latest_survival_time + 2, MAX_FOLLOW_UP, by=2)
    }else{
      seq(cur_visit_time, MAX_FOLLOW_UP, by=2)
    }
    res$biennial = getConsequences(SURV_CACHE_TIMES, SURV_CACHE_FULL, MAX_FOLLOW_UP,
                                   biennial_schedule, latest_survival_time)
    
    getPRIASSchedule(latest_survival_time, 
                     patient_data$psa[!is.na(patient_data$psa)],
                     patient_data$year_visit[!is.na(patient_data$psa)],
                     visit_schedule,
                     apply(PSA_CACHE_FULL, MARGIN = 1, sample, size=1))
    #now lets make the PRIAS schedule
    #But to account for variation in future PSA, we try it 50 times
    prias_conseq = lapply(1:50, function(x){
      prias_schedule = getPRIASSchedule(latest_survival_time, 
                                        patient_data$psa[!is.na(patient_data$psa)],
                                        patient_data$year_visit[!is.na(patient_data$psa)],
                                        visit_schedule,
                                        apply(PSA_CACHE_FULL, MARGIN = 1, sample, size=1))
      return(getConsequences(SURV_CACHE_TIMES, SURV_CACHE_FULL,MAX_FOLLOW_UP,
                             prias_schedule, latest_survival_time))
    })
    
    mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    prias_expected_delays = sapply(prias_conseq, FUN = "[[", "expected_delay")
    mode_expected_delay = mode(prias_expected_delays)
    mode_prias_conseq = which(sapply(prias_conseq, FUN = function(x){
      return(x$expected_delay == mode_expected_delay)
    })==T)[1]
    res$prias = prias_conseq[[mode_prias_conseq]]
  }
  
  return(list(schedules=res, horizon_surv_prob=horizon_surv_prob))
}

getPRIASSchedule = function(latest_survival_time, obs_psa, obs_psa_times, 
                            visit_schedule, future_log2psaplus1){
  
  fixed_schedule = c(1,4,7,10,15)
  
  proposed_biopsy_times = c()
  latest_biopsy_time = latest_survival_time
  
  for(i in 1:length(visit_schedule)){
    log2psa = log(c(obs_psa, pmax(0.1, 2^(future_log2psaplus1[1:i])-1)), base = 2)
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
    if(mean(surv_schedule_temp[i,], na.rm = T) <= surv_threshold & 
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
  
  #proposed biopsy times always has a final biopsy at year 10
  return(list(expected_delay = expected_delay,
              total_biopsies = length(orig_biopsy_times),
              biopsy_times=orig_biopsy_times,
              proposed_biopsy_times = proposed_biopsy_times))
}

