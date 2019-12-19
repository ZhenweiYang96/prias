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

times = allBiopsySchedules(cur_visit_time = 1, min_biopsy_gap = 1, horizon = 10)
print(all(!sapply(times, FUN = function(x){any(diff(x)<1)})))
