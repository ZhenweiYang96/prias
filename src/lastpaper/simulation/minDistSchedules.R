args = commandArgs(trailingOnly = F)

sim_seed = as.numeric(args[1])

library(JMbayes)
library(splines)
load(paste0("Rdata/lastpaper/sims/sim_seed_", sim_seed, ".Rdata"))
source("src/lastpaper/prediction_only_psa.R")
source("src/lastpaper/minDistThreshold.R")
source("src/lastpaper/scheduleCreatorCacheBased.R")

timesPerSubject = max(sim_res$testData$testDs$visitNumber)

decision_list = vector("list", length=timesPerSubject)

set.seed(sim_seed)
for(nb1_to_delay in seq(0, 2, 0.1)){
  nbCol = paste0("nb_", nb1_to_delay)
  delayCol = paste0("delay_", nb1_to_delay)
  
  sim_res$testData$testDs.id[, nbCol] = 0
  sim_res$testData$testDs.id[, delayCol] = NA
  
  for(testId in sim_res$testData$testDs.id$P_ID){
    
    print(paste("Doing for patient:", testId))
    curVisitNr = 5
    patient_data = sim_res$testData$testDs[sim_res$testData$testDs$P_ID == testId,]
    patient_data$gleason_sum[1] = 6
    
    idFilter = sim_res$testData$testDs.id$P_ID == testId
    progression_time = sim_res$testData$testDs.id$progression_time[idFilter]
    
    latest_survival_time = 0
    
    while(curVisitNr <= timesPerSubject){
      cur_visit_time = patient_data$year_visit[curVisitNr]
      if(cur_visit_time - latest_survival_time >= 1){
        decision = minDistScheduleDecision(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], 
                                           cur_visit_time = cur_visit_time,
                                           latest_survival_time = latest_survival_time,
                                           nb1_to_delay=nb1_to_delay,
                                           weight_by_horizon_risk = F)
        
        decision_list[[curVisitNr]] = decision
        if(decision$result==T){
          patient_data$gleason_sum[curVisitNr] = 6
          sim_res$testData$testDs.id[idFilter, nbCol] = sim_res$testData$testDs.id[idFilter, nbCol] + 1
          sim_res$testData$testDs.id[idFilter, delayCol] = cur_visit_time - progression_time
          latest_survival_time = cur_visit_time
          if(sim_res$testData$testDs.id[idFilter, delayCol]>=0){
            break;
          }
        }
      }
      
      curVisitNr = curVisitNr + 1
    }
  }
  
  save.image("Rdata/lastpaper/sim_seed_2019_sched_wtbyhorizonriskFALSE.Rdata")
}

decisions_pt25 = lapply(9:20, FUN = function(curVisitNr){
  decision = minDistScheduleDecision(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], 
                                     cur_visit_time = patient_data$year_visit[curVisitNr],
                                     latest_survival_time = 1,
                                     nb1_to_delay=0.25,
                                     weight_by_horizon_risk = F)
  
  return(decision)
})

patient_data = sim_res$testData$testDs[sim_res$testData$testDs$P_ID == 905,]
patient_data$gleason_sum[1] = 6

decisions_pt5 = lapply(9:20, FUN = function(curVisitNr){
  decision = minDistScheduleDecision(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], 
                                     cur_visit_time = patient_data$year_visit[curVisitNr],
                                     latest_survival_time = 1,
                                     nb1_to_delay=0.5,
                                     weight_by_horizon_risk = F)
  
  return(decision)
})

decisions_1 = lapply(9:20, FUN = function(curVisitNr){
  decision = minDistScheduleDecision(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], 
                                     cur_visit_time = patient_data$year_visit[curVisitNr],
                                     latest_survival_time = 1,
                                     nb1_to_delay=1,
                                     weight_by_horizon_risk = F)
  
  return(decision)
})

decisions_1pt5 = lapply(9:20, FUN = function(curVisitNr){
  decision = minDistScheduleDecision(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], 
                                     cur_visit_time = patient_data$year_visit[curVisitNr],
                                     latest_survival_time = 1,
                                     nb1_to_delay=1.5,
                                     weight_by_horizon_risk = F)
  
  return(decision)
})

sapply(1:length(decisions_1), function(i){
  
  x = decisions_1pt5[[i]]
  temp = seq(1,length(x$risk_thresholds), by = 100)
  
  pp=ggplot() + geom_line(aes(y=x$dist, x=x$risk_thresholds)) + 
    #geom_label(aes(y=x$dist[temp],x=x$risk_thresholds[temp], 
    #               label=paste0(round(x$expected_delays[temp],1),"--", x$total_biopsies[temp]))) +
  ggtitle(paste(x$latest_survival_time, "----",x$cur_visit_time)) + xlim(0,1)  
  
  y = decisions_pt25[[i]]
  pp2=ggplot() + geom_line(aes(y=y$dist, x=y$risk_thresholds)) + 
    #geom_label(aes(y=y$dist[temp],x=y$risk_thresholds[temp], 
    #               label=paste0(round(y$expected_delays[temp],1),"--", y$total_biopsies[temp]))) +
    ggtitle(paste(y$latest_survival_time, "----",y$cur_visit_time)) + xlim(0,1)  
  print(ggpubr::ggarrange(pp, pp2, ncol = 2, nrow = 1))
})

