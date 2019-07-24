timesPerSubject = max(sim_res$testData$testDs$visitNumber)

sim_res$testData$testDs$gleason_sum = NA
sim_res$testData$testDs.id$nb_mindist = 0
sim_res$testData$testDs.id$delay_mindist = NA
sim_res$testData$testDs$psa = 2^sim_res$testData$testDs$log2psaplus1 - 1
sim_res$trainingData$trainingDs$right_cens_time = sim_res$trainingData$trainingDs$progression_time
sim_res$trainingData$trainingDs$reclassification = sim_res$trainingData$trainingDs$progressed==1

PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP, 0.5))
#DRE check up time years
DRE_CHECK_UP_TIME = seq(0, MAX_FOLLOW_UP, 0.5)
BIOPSY_TEST_TIMES = DRE_CHECK_UP_TIME

roc_results_cache = vector("list", length(BIOPSY_TEST_TIMES))
names(roc_results_cache) = BIOPSY_TEST_TIMES
lapply(1:length(BIOPSY_TEST_TIMES), function(i){
  roc_results_cache[[i]] <<- vector("list", length(BIOPSY_TEST_TIMES))
  names(roc_results_cache[[i]]) <<- BIOPSY_TEST_TIMES
})

#Risk based schedule
# for(threshold in c(0.1, 0.05, 0.15)){
for(threshold in c(0.05, 0.15)){
  
  nbCol = paste0("nb_", threshold)
  delayCol = paste0("delay_", threshold)
  
  sim_res$testData$testDs.id[, nbCol] = 0
  sim_res$testData$testDs.id[, delayCol] = NA
  
  for(testId in sim_res$testData$testDs.id$P_ID){
    
    print(paste("Doing for patient:", testId))
    curVisitNr = 5
    patient_data = sim_res$testData$testDs[sim_res$testData$testDs$P_ID == testId,]
    patient_data$gleason_sum[1] = 6
    
    idFilter = sim_res$testData$testDs.id$P_ID == testId
    progressed = sim_res$testData$testDs.id$progressed[idFilter]
    progression_time = sim_res$testData$testDs.id$progression_time[idFilter]
    
    latest_survival_time = 0
    
    while(curVisitNr <= timesPerSubject){
      cur_visit_time = patient_data$year_visit[curVisitNr]
      
      if(cur_visit_time - latest_survival_time >= 1){
        set.seed(2019)
        gof = roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]]
        if(is.null(gof)){
          gof = goodness_of_fit(sim_res$mvJoint_psa_simDs, sim_res$trainingData$trainingDs,
                                T_start = latest_survival_time, T_horiz = cur_visit_time, horizon = MAX_FOLLOW_UP, M = 800)
          roc_results_cache[[as.character(latest_survival_time)]][[as.character(cur_visit_time)]] = gof
        }
        
        cumRiskMeans = rowMeans(1-getExpectedFutureOutcomes(sim_res$mvJoint_psa_simDs, patient_data[1:curVisitNr,], latest_survival_time, Inf,
                                                            survival_predict_times = c(cur_visit_time, gof$composite_horiz))$predicted_surv_prob)
        
        temp = gof$roc_results_composite[which(abs(gof$roc_results_composite$fpr - gof$roc_results$fpr[gof$roc_results$threshold==threshold])<= 1e-4 &
                                                 gof$roc_results_composite$tpr > gof$roc_results$tpr[gof$roc_results$threshold==threshold]),]
        if(nrow(temp)>0){
          temp = temp[order(temp$tpr, decreasing = T),]
          decision = cumRiskMeans[1] > temp$threshold1[1] &
            cumRiskMeans[2] > temp$threshold2[1]
        }else{
          decision = cumRiskMeans[1] >= threshold
        }
        
        if(decision==T){
          patient_data$gleason_sum[curVisitNr] = 6
          sim_res$testData$testDs.id[idFilter, nbCol] = sim_res$testData$testDs.id[idFilter, nbCol] + 1
          sim_res$testData$testDs.id[idFilter, delayCol] = cur_visit_time - progression_time
          latest_survival_time = cur_visit_time
          if(progressed & sim_res$testData$testDs.id[idFilter, delayCol]>=0){
            break;
          }
        }
      }
      
      curVisitNr = curVisitNr + 1
    }
  }
  
  testDs.id = sim_res$testData$testDs.id
  save(testDs.id, file="Rdata/lastpaper/sim_seed_2019_gof.Rdata")
}
