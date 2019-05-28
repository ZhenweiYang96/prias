source('src/clinical_gap3/constants.R')
source('src/clinical_gap3/prediction_only_psa.R')

compareSchedules = function(patient_data, cur_visit_time=NA, latest_survival_time=NA,
                            risk_thresholds = c(0.05, 0.1, 0.15), 
                            horizon = 10, 
                            M=500, CACHE_SIZE=200){
  #The various schedules we will compare are
  #biopsy every year
  #biopsy every 2 years
  #PRIAS
  #risk 5%, 10% and 15%
  
  if(is.na(cur_visit_time)){
    cur_visit_time = max(patient_data$year_visit)
    cur_visit_time = min(PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME >= cur_visit_time], horizon)
  }
  
  patient_data = patient_data[patient_data$year_visit <= cur_visit_time,]
  
  if(is.na(latest_survival_time)){
    latest_survival_time = max(patient_data$year_visit[!is.na(patient_data$gleason_sum)])
  }
  
  if(cur_visit_time < latest_survival_time){
    stop("Current visit time should be more than latest survival time")
  }else if((cur_visit_time - latest_survival_time)<MIN_BIOPSY_GAP){
    stop("No need for biopsy as current time should be more than latest survival by the minimum biopsy gap")
  }
  
  SURV_CACHE_TIMES = c(latest_survival_time, seq(latest_survival_time+1/365, horizon, length.out = CACHE_SIZE-1))
  pred_res=getExpectedFutureOutcomes(mvJoint_psa, patient_data, latest_survival_time, 
                                     survival_predict_times = SURV_CACHE_TIMES[-1],
                                     psa_predict_times = SURV_CACHE_TIMES, 
                                     psaDist = "Tdist", M = M, addRandomError = T)
  
  SURV_CACHE_PSA = pred_res$predicted_psa
  SURV_CACHE_MEAN = c(1,rowMeans(pred_res$predicted_surv_prob))
  #As many rows as the times at which surv is calculated
  SURV_CACHE_FULL = rbind(rep(1, M), pred_res$predicted_surv_prob)
  rm(pred_res)
  
  #First convert to risk
  COND_SURV_CACHE_MEAN = 1-SURV_CACHE_FULL
  #then scale the probabilities
  COND_SURV_CACHE_MEAN = COND_SURV_CACHE_MEAN / matrix(COND_SURV_CACHE_MEAN[CACHE_SIZE,], nrow=CACHE_SIZE, ncol=M, byrow = T)
  COND_SURV_CACHE_MEAN = apply(1-COND_SURV_CACHE_MEAN,1, mean)
  
  horizon_surv_prob = min(SURV_CACHE_MEAN)
  
  visit_schedule = PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME >= cur_visit_time & PSA_CHECK_UP_TIME<=horizon]
  visit_nearest_indices = sapply(visit_schedule, FUN = function(x){
    which.min(abs(x-SURV_CACHE_TIMES))
  })  
  surv_visit_schedule = SURV_CACHE_MEAN[visit_nearest_indices]
  SURV_CACHE_PSA = SURV_CACHE_PSA[visit_nearest_indices,]
  
  #Now we calculate schedule and consequences
  
  #added 3 for prias, annual and biennial schedule
  res = vector("list", length=length(risk_thresholds)+ 3)
  names(res) = c(paste("Risk:", risk_thresholds), "prias","annual", "biennial")
  for(i in 1:length(risk_thresholds)){
    threshold = risk_thresholds[i]
    risk_schedule = getRiskSchedule(SURV_CACHE_TIMES, SURV_CACHE_MEAN, threshold, horizon, 
                                    latest_survival_time, visit_schedule, surv_visit_schedule)
    res[[i]] = getConsequences(SURV_CACHE_TIMES, COND_SURV_CACHE_MEAN,
                               risk_schedule, latest_survival_time)
  }
  
  #first make the annual schedule
  annual_schedule = seq(cur_visit_time, horizon, by=1)
  if(max(annual_schedule)>=(horizon-0.5)){
    annual_schedule = annual_schedule[-length(annual_schedule)]
  }
  annual_schedule = c(annual_schedule, horizon)
  res$annual = getConsequences(SURV_CACHE_TIMES, COND_SURV_CACHE_MEAN,
                               annual_schedule, latest_survival_time)
  
  #now make the biennial schedule
  biennial_schedule = seq(cur_visit_time, horizon, by=2)
  if(max(biennial_schedule)>=(horizon-0.5)){
    biennial_schedule = biennial_schedule[-length(biennial_schedule)]
  }
  biennial_schedule = c(biennial_schedule, horizon)
  res$biennial = getConsequences(SURV_CACHE_TIMES, COND_SURV_CACHE_MEAN,
                                 biennial_schedule, latest_survival_time)
  
  #now lets make the PRIAS schedule
  #But to account for variation in future PSA, we try it 50 times
  prias_conseq = lapply(1:50, function(x){
    prias_schedule = getPRIASSchedule(latest_survival_time, 
                                      patient_data$psa[!is.na(patient_data$psa)],
                                      patient_data$year_visit[!is.na(patient_data$psa)],
                                      visit_schedule,
                                      apply(SURV_CACHE_PSA, MARGIN = 1, sample, size=1))
    if(max(prias_schedule)>=(horizon-0.5)){
      prias_schedule = prias_schedule[-length(prias_schedule)]
    }
    prias_schedule = c(prias_schedule, horizon)
    return(getConsequences(SURV_CACHE_TIMES, COND_SURV_CACHE_MEAN,
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

  return(list(schedules=res, horizon_surv_prob=horizon_surv_prob))
}

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
      if((visit_schedule[i] - latest_biopsy_time) >= MIN_BIOPSY_GAP){
        latest_biopsy_time = visit_schedule[i]
        proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
      }
      #else wait
    }else{
      if(visit_schedule[i] %in% fixed_schedule){
        if((visit_schedule[i] - latest_biopsy_time) >= MIN_BIOPSY_GAP){
          latest_biopsy_time = visit_schedule[i]
          proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
        }
        #else wait
      }
    }
  }
  
  return(proposed_biopsy_times)
}

getRiskSchedule = function(SURV_CACHE_TIMES, SURV_CACHE_MEAN,  
                           risk_threshold, horizon,
                           latest_survival_time,
                           visit_schedule, surv_visit_schedule){
  proposed_biopsy_times = c()
  latest_surv_prob = 1
  latest_biopsy_time = latest_survival_time
  for(i in 1:length(visit_schedule)){
    if((latest_surv_prob - surv_visit_schedule[i]) >= risk_threshold & 
       (visit_schedule[i]-latest_biopsy_time)>=MIN_BIOPSY_GAP){
      latest_biopsy_time = visit_schedule[i]
      proposed_biopsy_times = c(proposed_biopsy_times, visit_schedule[i])
      latest_surv_prob = surv_visit_schedule[i]
    }
  }
  
  #Make the last biopsy at year 10
  #and if the last biopsy was 9.5 or more, remove it and only keep the one at 10
  if(max(proposed_biopsy_times)>=(horizon-0.5)){
    proposed_biopsy_times = proposed_biopsy_times[-length(proposed_biopsy_times)]
  }
  proposed_biopsy_times = c(proposed_biopsy_times, horizon)
  
  return(proposed_biopsy_times)
}

getConsequences = function(SURV_TIMES, SURV_MEAN, 
                           proposed_biopsy_times, latest_survival_time){
  
  #Now we create intervals in which to integrate the conditional surv prob
  if(length(proposed_biopsy_times)==1){
    intervals = list(c(latest_survival_time,proposed_biopsy_times))
  }else{
    intervals = c(latest_survival_time,
                  rep(proposed_biopsy_times[-length(proposed_biopsy_times)],each=2),
                  proposed_biopsy_times[length(proposed_biopsy_times)])
    intervals = split(intervals, rep(1:(length(intervals)/2), each=2))
  }
  
  delay_prob = lapply(intervals, function(interval){
    lower_limit = interval[1]
    upper_limit = interval[2]
    
    total_marginal_prob = SURV_MEAN[which.min(abs(lower_limit-SURV_TIMES))] -
      SURV_MEAN[which.min(abs(upper_limit-SURV_TIMES))]
    
    wt_points=getGaussianQuadWeightsPoints(interval)
    
    expected_fail_time = lower_limit + sum(sapply(wt_points$points, FUN = function(x){
      prob_at_x_scaled = SURV_MEAN[which.min(abs(x-SURV_TIMES))] -
        SURV_MEAN[which.min(abs(upper_limit-SURV_TIMES))]
      prob_at_x_scaled = prob_at_x_scaled / total_marginal_prob
    }) * wt_points$weights)
    
    delay = upper_limit - expected_fail_time
    
    return(c(delay=delay, total_marginal_prob=total_marginal_prob))
  })
  
  return(list(delay_prob = delay_prob,
              expected_delay = sum(sapply(delay_prob, FUN = prod)),
              prop_total_biopsies = length(proposed_biopsy_times),
              proposed_biopsy_times=proposed_biopsy_times))
}