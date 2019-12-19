seeds = 2101:2103

for(seed in seeds){
  files = list.files("Rdata/lastpaper/simulation/results/",
                     pattern = as.character(seed), full.names = T)
  
  biopsyDfs = lapply(files, FUN = function(file){
    load(file)
    return(biopsyDf)
  })
  
  biopsyDf_summary_list = lapply(biopsyDfs, FUN = function(patient_df){
    schedules = unique(patient_df$schedule)
    
    nb = sapply(schedules, FUN = function(schedule){
      nrow(patient_df[patient_df$schedule==schedule,])
    })
    
    delay = sapply(schedules, FUN = function(schedule){
      tail(patient_df$biopsy_time[patient_df$schedule==schedule],1) - patient_df$progression_time[1]
    })
    
    ret_df = patient_df[rep(1, length(schedules)), c("seed", "P_ID", "progression_time", "progressed")]
    ret_df$schedule = schedules
    ret_df$nb = nb
    ret_df$delay = delay
    return(ret_df)
  })
  
  biopsyDf_summary = do.call('rbind', biopsyDf_summary_list)
  
  save(biopsyDf_summary, file=paste0("Rdata/lastpaper/simulation/combined_results/seed_", seed, ".Rdata"))
}
