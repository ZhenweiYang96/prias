createSchedules = function(patient_data, cur_visit_time, latest_survival_time,
                            risk_thresholds = c(0.05, 0.1, 0.15), 
                            min_biopsy_gap, SURV_CACHE_TIMES, PSA_CACHE_TIMES,
                            SURV_CACHE_FULL, PSA_CACHE_FULL){
  #The various schedules we will compare are
  #biopsy every year
  #biopsy every 2 years
  #PRIAS
  #risk 5%, 10% and 15%
  
  # if(cur_visit_time - latest_survival_time < min_biopsy_gap){
  #   cur_visit_time = latest_survival_time + min_biopsy_gap
  # }
  
  #Decision epochs in years
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP, 0.5))
  
  horizon_surv_prob = min(rowMeans(SURV_CACHE_FULL))
  
  visit_schedule = c(cur_visit_time, PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > cur_visit_time & PSA_CHECK_UP_TIME<=MAX_FOLLOW_UP])
  visit_nearest_indices = sapply(visit_schedule, FUN = function(x){
    which.min(abs(x-SURV_CACHE_TIMES))
  })
  surv_schedule = SURV_CACHE_FULL[visit_nearest_indices,]
  PSA_CACHE_FULL = PSA_CACHE_FULL[visit_nearest_indices,]
  
  #added 3 for prias, annual and biennial schedule
  res = vector("list", length=length(risk_thresholds)+ 3)
  names(res) = c(paste("Risk:", risk_thresholds), "prias","annual", "biennial")
  
  for(i in 1:length(risk_thresholds)){
    threshold = risk_thresholds[i]
    risk_schedule = getRiskSchedule(1-threshold, MAX_FOLLOW_UP, 
                                    latest_survival_time, visit_schedule, surv_schedule, min_biopsy_gap)
    res[[i]] = getConsequences(SURV_CACHE_TIMES, SURV_CACHE_FULL,
                               risk_schedule, latest_survival_time)
  }
  
  #first make the annual schedule
  annual_schedule = if(cur_visit_time - latest_survival_time < min_biopsy_gap){
    seq(latest_survival_time + min_biopsy_gap, MAX_FOLLOW_UP, by=1)
  }else{
    seq(cur_visit_time, MAX_FOLLOW_UP, by=1)
  }
  if(max(annual_schedule)>=(MAX_FOLLOW_UP-0.5)){
    annual_schedule = annual_schedule[-length(annual_schedule)]
  }
  annual_schedule = c(annual_schedule, MAX_FOLLOW_UP)
  res$annual = getConsequences(SURV_CACHE_TIMES, SURV_CACHE_FULL,
                               annual_schedule, latest_survival_time)
  
  #now make the biennial schedule
  biennial_schedule = if(cur_visit_time - latest_survival_time < 2){
    seq(latest_survival_time + 2, MAX_FOLLOW_UP, by=2)
  }else{
    seq(cur_visit_time, MAX_FOLLOW_UP, by=2)
  }
  if(max(biennial_schedule)>=(MAX_FOLLOW_UP-0.5)){
    biennial_schedule = biennial_schedule[-length(biennial_schedule)]
  }
  biennial_schedule = c(biennial_schedule, MAX_FOLLOW_UP)
  res$biennial = getConsequences(SURV_CACHE_TIMES, SURV_CACHE_FULL,
                                 biennial_schedule, latest_survival_time)
  
  #now lets make the PRIAS schedule
  #But to account for variation in future PSA, we try it 50 times
  prias_conseq = lapply(1:50, function(x){
    prias_schedule = getPRIASSchedule(latest_survival_time, 
                                      patient_data$psa[!is.na(patient_data$psa)],
                                      patient_data$year_visit[!is.na(patient_data$psa)],
                                      visit_schedule,
                                      apply(PSA_CACHE_FULL, MARGIN = 1, sample, size=1))
    if(max(prias_schedule)>=(MAX_FOLLOW_UP-0.5)){
      prias_schedule = prias_schedule[-length(prias_schedule)]
    }
    prias_schedule = c(prias_schedule, MAX_FOLLOW_UP)
    return(getConsequences(SURV_CACHE_TIMES, SURV_CACHE_FULL,
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
