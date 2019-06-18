compareSchedules = function(patient_data, cur_visit_time=NA, latest_survival_time=NA,
                            risk_thresholds = c(0.05, 0.1, 0.15), 
                            horizon = 10, weight_by_horizon_risk=T,
                            min_biopsy_gap = 1,
                            M=500, CACHE_SIZE=200){
  #The various schedules we will compare are
  #biopsy every year
  #biopsy every 2 years
  #PRIAS
  #risk 5%, 10% and 15%
  
  #Decision Epochs in years
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
  #DRE check up time years
  DRE_CHECK_UP_TIME = seq(0, horizon, 0.5)
  
  BIOPSY_TEST_TIMES = DRE_CHECK_UP_TIME
  
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
  }else if((cur_visit_time - latest_survival_time)<min_biopsy_gap){
    stop("No need for biopsy as current time should be more than latest survival by the minimum biopsy gap")
  }
  
  SURV_CACHE_TIMES = c(latest_survival_time, seq(latest_survival_time+1/365, horizon, length.out = CACHE_SIZE-1))
  pred_res=getExpectedFutureOutcomes(mvJoint_psa_time_scaled, patient_data, latest_survival_time, 
                                     survival_predict_times = SURV_CACHE_TIMES[-1],
                                     psa_predict_times = SURV_CACHE_TIMES, 
                                     psaDist = "Tdist", M = M, addRandomError = T)
  
  PSA_CACHE_FULL = pred_res$predicted_psa
  SURV_CACHE_MEAN = c(1,rowMeans(pred_res$predicted_surv_prob))
  #As many rows as the times at which surv is calculated
  SURV_CACHE_FULL = rbind(rep(1, M), pred_res$predicted_surv_prob)
  rm(pred_res)
  
  #First convert to conditional risk (condition on T<horizon)
  COND_RISK_CACHE_FULL = t(apply(1-SURV_CACHE_FULL, 1, FUN = function(x){
    x / (1-SURV_CACHE_FULL[CACHE_SIZE,])
  }))
  COND_SURV_CACHE_FULL = 1-COND_RISK_CACHE_FULL
  #See I could do this 1- thing in one step, but this is for readability
  rm(COND_RISK_CACHE_FULL)
  
  if(weight_by_horizon_risk==T){
    COND_SURV_CACHE_FULL = SURV_CACHE_FULL
  }
  
  horizon_surv_prob = min(SURV_CACHE_MEAN)
  
  visit_schedule = PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME >= cur_visit_time & PSA_CHECK_UP_TIME<=horizon]
  visit_nearest_indices = sapply(visit_schedule, FUN = function(x){
    which.min(abs(x-SURV_CACHE_TIMES))
  })  
  surv_schedule = SURV_CACHE_FULL[visit_nearest_indices,]
  PSA_CACHE_FULL = PSA_CACHE_FULL[visit_nearest_indices,]
  
  #Now we calculate schedule and consequences
  
  #added 3 for prias, annual and biennial schedule
  res = vector("list", length=length(risk_thresholds)+ 3)
  names(res) = c(paste("Risk:", risk_thresholds), "prias","annual", "biennial")
  
  for(i in 1:length(risk_thresholds)){
    threshold = risk_thresholds[i]
    risk_schedule = getRiskSchedule(1-threshold, horizon, 
                                    latest_survival_time, visit_schedule, surv_schedule, min_biopsy_gap)
    res[[i]] = getConsequences(SURV_CACHE_TIMES, COND_SURV_CACHE_FULL,
                               risk_schedule, latest_survival_time)
  }
  
  #first make the annual schedule
  annual_schedule = seq(cur_visit_time, horizon, by=1)
  if(max(annual_schedule)>=(horizon-0.5)){
    annual_schedule = annual_schedule[-length(annual_schedule)]
  }
  annual_schedule = c(annual_schedule, horizon)
  res$annual = getConsequences(SURV_CACHE_TIMES, COND_SURV_CACHE_FULL,
                               annual_schedule, latest_survival_time)
  
  #now make the biennial schedule
  biennial_schedule = seq(cur_visit_time, horizon, by=2)
  if(max(biennial_schedule)>=(horizon-0.5)){
    biennial_schedule = biennial_schedule[-length(biennial_schedule)]
  }
  biennial_schedule = c(biennial_schedule, horizon)
  res$biennial = getConsequences(SURV_CACHE_TIMES, COND_SURV_CACHE_FULL,
                                 biennial_schedule, latest_survival_time)
  
  #now lets make the PRIAS schedule
  #But to account for variation in future PSA, we try it 50 times
  prias_conseq = lapply(1:50, function(x){
    prias_schedule = getPRIASSchedule(latest_survival_time, 
                                      patient_data$psa[!is.na(patient_data$psa)],
                                      patient_data$year_visit[!is.na(patient_data$psa)],
                                      visit_schedule,
                                      apply(PSA_CACHE_FULL, MARGIN = 1, sample, size=1),
                                      min_biopsy_gap)
    if(max(prias_schedule)>=(horizon-0.5)){
      prias_schedule = prias_schedule[-length(prias_schedule)]
    }
    prias_schedule = c(prias_schedule, horizon)
    return(getConsequences(SURV_CACHE_TIMES, COND_SURV_CACHE_FULL,
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
source("src/clinical_gap3/schedules.R")