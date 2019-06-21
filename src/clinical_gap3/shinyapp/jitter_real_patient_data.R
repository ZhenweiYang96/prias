demo_pat_list = list(pat1_data, pat2_data, pat3_data, pat4_data)
demo_pat_list = lapply(1:length(demo_pat_list), function(i){
  set.seed(2019 + i)
  x = demo_pat_list[[i]][,c("P_ID", "age", "start_date", "year_visit", "psa", "gleason_sum")]
  x$age = x$age + rnorm(1, 0, 2)
  x$psa = x$psa + rnorm(nrow(x), 0, 1)
  gleason_sums = x$gleason_sum[!is.na(x$gleason_sum)]
  gleason_sums = sapply(gleason_sums, FUN = function(gs){
    if(gs<=6){
      gs = 6
    }else if(gs>=7){
      gs = 7
    }
  })
  x$gleason_sum[!is.na(x$gleason_sum)] = gleason_sums
  min_diff_year_visit = min(diff(x$year_visit))
  x$year_visit[-1] = x$year_visit[-1] + runif(n=length(x$year_visit[-1]),
                                              min = 1/365, 
                                              max = min_diff_year_visit)
  return(x)
})