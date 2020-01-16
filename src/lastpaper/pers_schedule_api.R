#Only works for one patient at a time
personalizedSchedule.mvJMbayes <- function (object, newdata, idVar = "id", last_test_time = NULL,
                                            gap = 0, horizon, detection_delay_limit=Inf, total_tests_limit=Inf,
                                            seed = 1L, M = 200L, cache_size=1000, ...) {
  
  if(length(unique(newdata[[idVar]])) > 1){
    stop("Personalized schedule can be created for only one patient at a time.\n")
  }
  
  timeVar <- object$model_info$timeVar
  
  if(is.null(last_test_time)){
    last_test_time <- max(newdata[[timeVar]], na.rm = T)
  }
  
  cur_visit_time = max(newdata[[timeVar]], na.rm = T)
  
  #Step 1: Create a CACHE of predicted survival probabilities to speed up operations
  SURV_CACHE_TIMES <- seq(from = last_test_time, to = horizon, length.out = cache_size)
  
  surv_pred <- survfitJM(object=object, newdata=newdata, idVar=idVar, survTimes = SURV_CACHE_TIMES[-1],
                         last.time = last_test_time, seed=seed, M=M, ...)
  
  SURV_CACHE_FULL <- rbind(rep(1, M), 
                           do.call('cbind', lapply(surv_pred$full.results, "[[", 1)))
  
  ##############################################
  #Definition of the function for fixed schedule
  ##############################################
  fixedSchedule <- function(surv_threshold){
    proposed_test_times <- c()
    previous_test_time <- last_test_time
    surv_schedule_temp <- SURV_CACHE_FULL
    
    for(i in 1:length(SURV_CACHE_TIMES)){
      if(SURV_CACHE_TIMES[i] - previous_test_time>=gap){
        if(mean(surv_schedule_temp[i,], na.rm = T) <= surv_threshold){
          previous_test_time <- SURV_CACHE_TIMES[i]
          proposed_test_times <- c(proposed_test_times, previous_test_time)
          surv_schedule_temp <- apply(SURV_CACHE_FULL, 2, FUN = function(x){x/x[i]})
        }
      }
    }
    
    return(proposed_test_times)
  }
  
  wk = JMbayes:::gaussKronrod()$wk
  sk = JMbayes:::gaussKronrod()$sk
  getGaussianQuadWeightsPoints = function(lower_upper_limit){
    p1 = 0.5 * sum(lower_upper_limit)
    p2 = 0.5 * diff(lower_upper_limit)
    
    return(list(weights = p2*wk, points = p2 * sk + p1))
  }
  
  ###############################################
  # Consequences
  ###############################################
  consequences = function(proposed_test_times){
    
    practical_test_times = proposed_test_times
    #Make the last test at year 10
    if(horizon - tail(practical_test_times,1) < gap){
      practical_test_times = practical_test_times[-length(practical_test_times)]
    }
    practical_test_times = c(practical_test_times, horizon)
    
    #Now we create intervals in which to integrate the conditional surv prob
    if(length(practical_test_times)==1){
      test_intervals = list(c(last_test_time, practical_test_times))
    }else{
      test_intervals = c(last_test_time,
                         rep(practical_test_times[-length(practical_test_times)],each=2),
                         practical_test_times[length(practical_test_times)])
      test_intervals = split(test_intervals, rep(1:(length(test_intervals)/2), each=2))
    }
    
    interval_res = vector("list", length(test_intervals))
    for(j in 1:length(test_intervals)){
      wt_points=getGaussianQuadWeightsPoints(test_intervals[[j]])
      lower_limit = test_intervals[[j]][1]
      upper_limit = test_intervals[[j]][2]
      
      lower_limit_nearest_index = which.min(abs(lower_limit-SURV_CACHE_TIMES))
      upper_limit_nearest_index = which.min(abs(upper_limit-SURV_CACHE_TIMES))
      interval_res[[j]]$cum_risk_interval = SURV_CACHE_FULL[lower_limit_nearest_index,] - 
        SURV_CACHE_FULL[upper_limit_nearest_index,]
      
      cond_expected_fail_time = sapply(1:length(wt_points$points), function(i){
        cum_surv_at_points = SURV_CACHE_FULL[which.min(abs(wt_points$points[i]-SURV_CACHE_TIMES)),] - SURV_CACHE_FULL[upper_limit_nearest_index,]
        scaled_cum_surv_at_points = cum_surv_at_points/interval_res[[j]]$cum_risk_interval
        
        return(wt_points$weights[i] * scaled_cum_surv_at_points)
      })
      
      interval_res[[j]]$delay = upper_limit - (lower_limit + apply(cond_expected_fail_time, 1, sum))
    }
    
    expected_delay = mean(apply(sapply(interval_res, FUN = function(x){
      x$delay * x$cum_risk_interval
    }),1, sum, na.rm=T),na.rm=T)
    
    #proposed test times always has a final test at year 10
    return(list(expected_delay = expected_delay,
                proposed_test_times=proposed_test_times,
                practical_test_times = practical_test_times))
  }
  
  #Now we create many fixed risk based schedules
  risk_thresholds <- seq(0, 1, length.out = 100)
  all_schedules <- vector("list", length=length(risk_thresholds))
  names(all_schedules) <- c(paste("Risk:", risk_thresholds))
  
  for(i in 1:length(risk_thresholds)){
    threshold <- risk_thresholds[i]
    risk_schedule <- fixedSchedule(1-threshold)
    if(is.null(risk_schedule)){
      risk_schedule <- horizon
    }
    
    all_schedules[[i]] = consequences(risk_schedule)
    all_schedules[[i]]$threshold = threshold
    all_schedules[[i]]$euclidean_distance = sqrt(all_schedules[[i]]$expected_delay^2 + (length(all_schedules[[i]]$practical_test_times)-1)^2)
  }
  
  expected_delays = sapply(all_schedules, "[[", "expected_delay")
  total_tests = sapply(lapply(all_schedules, "[[", "practical_test_times"), length)
  dist = sapply(all_schedules, "[[", "euclidean_distance")
  
  ############
  # this is the code for applying constraints before final optimization
  ###########
  flag = 0
  old_total_test_limit = total_tests_limit
  while(sum(total_tests<=total_tests_limit)==0){
    flag = 1
    total_tests_limit = total_tests_limit + 1
  }
  
  if(flag==1){
    print(paste0("Could not find any schedule with a total_tests_limit=", old_total_test_limit))
    print(paste0("Using new total_tests_limit=",total_tests_limit))
  }
  
  flag = 0
  old_detection_delay_limit = detection_delay_limit
  delta_delay = diff(range(expected_delays))/100
  while(sum(expected_delays<=detection_delay_limit) == 0){
    flag = 1
    detection_delay_limit = detection_delay_limit + delta_delay
  }
  
  if(flag==1){
    print(paste0("Could not find any schedule with a detection_delay_limit=",old_detection_delay_limit))
    print(paste0("Using new detection_delay_limit=",detection_delay_limit))
  }
  
  filter = expected_delays<=detection_delay_limit & total_tests<=total_tests_limit
  filtered_schedules = all_schedules[filter]
  filtered_dist = dist[filter]
  
  #There are many, but we select only one for now
  optimal_threshold_index = which.min(filtered_dist)[1]
  ret = filtered_schedules[[optimal_threshold_index]]
  ret$proposed_test_times = NULL
  names(ret) = c("expected_detection_delay", "test_schedule", "risk_threshold", "euclidean_distance")
  
  full_data = list('selected_schedule'=ret, 'all_schedules'=all_schedules)
  
  class(full_data) <- "personalizedSchedule.mvJMbayes"
  return(full_data)
}